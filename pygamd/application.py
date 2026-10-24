'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of PYGAMD do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with PYGAMD are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu, 
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''

# CuPy-native run loop driving the RawKernel forces / integrators.  No numba.
# The global zero-force step (non-fused path) is a CuPy array fill on the CuInfo
# buffers.  Optional particle sorting (dynamics(sort=True)) is supported via
# pygamd.tinker, imported lazily only on that path:
#   - the current sorter is a deterministic CELL-KEY sort (full Hilbert traversal
#     is partial/deferred);
#   - it refreshes tag/rtag and the bond/angle/dihedral index views through
#     CuInfo (tag-based topology is sort-invariant) and forces an nlist rebuild;
#   - sort=False remains the default and never imports tinker.
from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters (no-op unless PYGAMD_ENGINE_METRICS=1)
from pygamd import engine_timing as _t      # engine timing opt-in timing (no-op unless PYGAMD_ENGINE_TIMING=1)
from pygamd import launch_trace as _lt  # launch trace opt-in launch trace (no-op unless PYGAMD_LAUNCH_TRACE=1)
from pygamd.cuda_graph import StableCudaGraph
import numpy as np
import time
import math
import os


class dynamics:
	def __init__(self, info, dt, sort=False):
		self.forces = []
		self.integrations = []
		self.tinkers = []
		self.dumps = []

		self.info=info
		# CuPy-native main line: the snapshot provides a CuInfo device-data layer
		# (info.cu) that RawPair / RawNVT consume directly.  No CUDA driver work is
		# done here at construction; the global zero-force step is a CuPy fill.
		self.is_cuinfo = bool(getattr(info, "cu", None))
		self.info.dt = dt
		self.sort = sort
		# Eager-graph primary controller.  The default high-level policy is
		# PYGAMD_EXECUTION_POLICY=eager_graph; set eager_debug for the legacy pure-eager
		# debug route or graph_strict for CI/coverage enforcement.
		self._graph = StableCudaGraph()
		self.execution_policy = getattr(self._graph, "execution_policy", "eager_graph")
		# Application-loop bookkeeping.  Public run() now dispatches to an
		# explicit eager-graph primary loop or the demoted legacy eager_debug loop.
		# These fields are surfaced by engine diagnostics so route tests can
		# assert the application loop directly, not just the environment policy.
		self._application_loop_kind = None
		self._application_loop_label = None
		self._application_loop_count = 0
		self._application_loop_eager_graph_count = 0
		self._application_loop_eager_debug_count = 0
		self.start_time = time.perf_counter()
		self.init_timestep = self.info.timestep
		self.last = self.init_timestep
		# Default GALAMOST-compatible TPS target scheduling keeps the legacy
		# adaptive "about 20 seconds per report" behavior.  Benchmark drivers may
		# set PYGAMD_TPS_PERIOD_STEPS to force fixed step windows, while preserving
		# GALAMOST's aggregation rule that excludes the first window from Average TPS.
		self.tps_period_steps = int(os.environ.get("PYGAMD_TPS_PERIOD_STEPS", "0") or "0")
		if self.tps_period_steps < 0:
			self.tps_period_steps = 0
		self.tps_period_mode = "fixed_steps" if self.tps_period_steps > 0 else "adaptive_20s"
		self.target = self.init_timestep + (self.tps_period_steps if self.tps_period_steps > 0 else 200)
		self.begin = self.init_timestep 
		self.end = self.init_timestep
		# TPS accounting follows GALAMOST Application::TPScompute(): the
		# first TPS window is treated as a warm-up/sample-alignment window and is
		# not included in Average TPS.  The cumulative counters intentionally span
		# segmented runs such as app.run(warmup); app.run(timed), matching
		# GALAMOST's application-level rule instead of an external shell stopwatch.
		self.average_tps = 0.0
		self.count_tps = 0
		self.first_tps = True
		self.tps_samples = []
		self.average_tps_value = None
		self.tps_window_status = "INSUFFICIENT_TPS_WINDOW"
		self.block_size = 64
		if dt>0.0:
			self.sort_period = np.int32(2.0/dt)
			if self.sort_period<100:
				self.sort_period = 100
			elif self.sort_period > 2000:
				self.sort_period = 2000
		else:
			self.sort_period = np.int32(1000)
		if self.sort:
			from pygamd import tinker   # lazy: CuPy-native cell-key sorter, only when sorting
			self.psort = tinker.psort(info, self.sort_period)

	def run(self, ntimesteps):
		"""Run dynamics using the resolved high-level execution policy.

		The graph-aware loop the normal application loop.  The legacy
		pure-eager implementation remains available only through the explicit
		eager_debug policy and lives in _run_eager_debug() so future cleanups can
		remove/demote pure-eager code without touching the graph-primary path.
		"""
		if self.execution_policy == "eager_debug":
			return self._run_eager_debug(ntimesteps)
		return self._run_eager_graph_primary(ntimesteps)

	def _select_application_loop(self, kind, label):
		self._application_loop_kind = str(kind)
		self._application_loop_label = str(label)
		self._application_loop_count += 1
		if kind == "eager_debug":
			self._application_loop_eager_debug_count += 1
		else:
			self._application_loop_eager_graph_count += 1
		# Mirror application-loop diagnostics onto the graph controller because
		# external collectors already pull their route metadata from it.
		# Compatibility aliases are kept while kind/label/count fields remain canonical, so older reports and newer route tests read the same state.
		graph_aware = (self._application_loop_kind != "eager_debug")
		setattr(self._graph, "application_loop_kind", self._application_loop_kind)
		setattr(self._graph, "application_loop_label", self._application_loop_label)
		setattr(self._graph, "application_loop_count", int(self._application_loop_count))
		setattr(self._graph, "application_loop_eager_graph_count", int(self._application_loop_eager_graph_count))
		setattr(self._graph, "application_loop_eager_debug_count", int(self._application_loop_eager_debug_count))
		setattr(self._graph, "application_loop", self._application_loop_kind)
		setattr(self._graph, "application_loop_graph_aware", bool(graph_aware))
		if _sp.ENABLED:
			_sp.set_meta("pygamd_application_loop_kind", self._application_loop_kind)
			_sp.set_meta("pygamd_application_loop_label", self._application_loop_label)
			_sp.set_meta("pygamd_application_loop_count", int(self._application_loop_count))
			_sp.set_meta("pygamd_application_loop_eager_graph_count", int(self._application_loop_eager_graph_count))
			_sp.set_meta("pygamd_application_loop_eager_debug_count", int(self._application_loop_eager_debug_count))

	def _run_begin_state(self, ntimesteps, loop_kind, loop_label):
		self._select_application_loop(loop_kind, loop_label)
		self.end += ntimesteps
		print("info : --- start to run")
		print("info : from ",self.begin," to ",self.end," timestep")

		# engine metrics: coarse run-level wall time of the WHOLE measured region of
		# app.run() (register + begin sort + begin force compute + dumps + the main
		# step loop), opt-in.  Started here -- before register -- so it covers the
		# same work the per-category timing accumulates, keeping wall_timed_region_ms
		# consistent with total_profiled_ms.  (Key name kept for the collector.)
		_metrics_loop_t0 = time.perf_counter() if _sp.ENABLED else None

		self.register(self.begin)
		if self.sort:
			self.psort.compute(self.begin)
		# Fused zero-force path: fuse global zero_force_vp into a single RawKernel pair force
		# when it is safe (see _fuse_target); otherwise keep the original
		# global-zero-then-additive behavior unchanged.
		#
		# Bug fix (segmented runs): the begin-timestep force compute must zero the
		# global force / virial_potential buffers on the ADDITIVE (non-fused) path,
		# exactly like the main loop does.  Without this, a second app.run() (e.g.
		# app.run(warmup); app.run(timed)) starts its begin force compute on top of
		# the previous run's force, polluting bonded / DPD / any multi-force result.
		# The fused single-RawPair path overwrites force in its kernel, so it must
		# NOT be zeroed here (preserves the fast-path).
		fused = self._apply_fuse_zero_force()
		if not fused:
			self._zero_force_vp()
		for fi in self.forces:
			fi.compute(self.begin)

		for di in self.dumps:
			di.compute(self.begin)
		return _metrics_loop_t0

	def _run_finish_state(self, _metrics_loop_t0):
		if _sp.ENABLED and _metrics_loop_t0 is not None:
			_sp.add_time("python_step_loop_wall_s", time.perf_counter() - _metrics_loop_t0)
		self.begin = self.end
		prof.finalize()
		self._graph.finalize()
		_sp.record_graph(self)   # engine metrics: capture CUDA-graph status/counts (opt-in)

	def _run_step_prefix(self, ts):
		if _sp.ENABLED:
			_sp.inc("step_loop_iters")
		if _t.ENABLED:
			_t.tick_step()
		if _lt.LAUNCH:
			_lt.tick_step()
		prof.begin("step_total")
		prof.count("steps")
		prof.begin("register")
		self.register(ts)
		if _lt.LAUNCH: _lt.action("register_calls")
		prof.end("register")

	def _run_integrator_firststep(self, ts):
		prof.begin("integration_firststep_total")
		for it in self.integrations:
			if _lt.LAUNCH: _lt.action("integrator_firststep_calls")
			it.data.firststep(ts)
		prof.end("integration_firststep_total")

	def _run_eager_force_integrator_tail(self, ts, fused, record_graph_fallback=False):
		"""Run the legacy eager force + integrator-tail segment.

		When called from the graph-primary loop, record_graph_fallback=True records a
		whole-step fallback for graph coverage.  The explicit eager_debug loop calls
		this same kernel substrate with record_graph_fallback=False, so pure eager is
		not misreported as graph fallback diagnostics.
		"""
		if record_graph_fallback:
			self._graph.record_eager_step()
		if not fused:
			prof.count("nonfused_zero_steps")
			prof.count("zero_force_vp_calls")
			prof.begin("zero_force_vp")
			self._zero_force_vp()
			prof.end("zero_force_vp")
		else:
			prof.count("fused_zero_steps")
			prof.count("zero_force_vp_skipped_fused")

		prof.begin("force_compute_total")
		for fi in self.forces:
			if _lt.LAUNCH: _lt.action("force_compute_calls")
			fi.compute(ts)
		prof.end("force_compute_total")

		prof.begin("integration_secondstep_total")
		for ii in self.integrations:
			if _lt.LAUNCH: _lt.action("integrator_secondstep_calls")
			ii.data.secondstep(ts)
		prof.end("integration_secondstep_total")

	def _run_step_suffix(self, ts):
		prof.begin("dump_compute_total")
		for di in self.dumps:
			if _lt.LAUNCH: _lt.action("dump_calls")
			di.compute(ts)
		prof.end("dump_compute_total")

		prof.begin("tinker_compute_total")
		for ti in self.tinkers:
			ti.data.compute(ts)
		prof.end("tinker_compute_total")
		if self.sort:
			prof.begin("sort_total")
			self.psort.compute(ts)
			prof.end("sort_total")
		self.tps_compute(ts)
		prof.end("step_total")
		prof.tick_step()

	def _run_eager_graph_primary(self, ntimesteps):
		_metrics_loop_t0 = self._run_begin_state(
			ntimesteps,
			"eager_graph_primary",
			"graph-aware primary application loop",
		)
		for ts in range(self.begin+1, self.end+1):
			self._run_step_prefix(ts)

			# Three-layer CUDA-graph stable segment.  FULL graph contains its OWN
			# firststep, so it must be decided BEFORE the eager firststep (otherwise
			# firststep would run twice).
			#   (1) FULL : firststep + force segment + integrator tail.
			#   (2) TAIL : eager firststep + eager nlist + graph force/tail.
			#   (3) EAGER: controlled fallback, counted for graph coverage.
			ran_graph = False
			if self._graph.enabled and self._graph.allow_full and self._graph.try_run_full(self, ts):
				ran_graph = True
				if _sp.ENABLED: _sp.inc("cuda_graph_full_step_count")
				prof.count("cuda_graph_full_step")

			if not ran_graph:
				self._run_integrator_firststep(ts)
				fused = self._apply_fuse_zero_force()

				# TAIL graph: run all graph-relevant nlists eagerly (authoritative) then
				# graph the force+integrator tail through StableCudaGraph.
				if self._graph.enabled and self._graph.prepare_tail_nlists(self, ts):
					if self._graph.try_run_tail(self, ts):
						ran_graph = True
						if _sp.ENABLED: _sp.inc("cuda_graph_tail_step_count")
						prof.count("cuda_graph_tail_step")

				if not ran_graph:
					self._run_eager_force_integrator_tail(ts, fused, record_graph_fallback=True)

			self._run_step_suffix(ts)
		self._run_finish_state(_metrics_loop_t0)

	def _run_eager_debug(self, ntimesteps):
		"""Legacy pure-eager debug/reference loop.

		This path intentionally bypasses full/tail graph attempts and does not record
		whole-step eager execution as graph fallback.  It keeps the same eager kernel
		substrate used by graph capture, but is no longer the public primary loop.
		"""
		_metrics_loop_t0 = self._run_begin_state(
			ntimesteps,
			"eager_debug",
			"legacy pure-eager debug application loop",
		)
		for ts in range(self.begin+1, self.end+1):
			self._run_step_prefix(ts)
			self._run_integrator_firststep(ts)
			fused = self._apply_fuse_zero_force()
			self._run_eager_force_integrator_tail(ts, fused, record_graph_fallback=False)
			self._run_step_suffix(ts)
		self._run_finish_state(_metrics_loop_t0)

	def add(self, app):
		if app.name=="force":
			self.forces.append(app)
		elif app.name=="integration":
			self.integrations.append(app)
		elif app.name=="tinker":
			self.tinkers.append(app)
		elif app.name=="dump":
			self.dumps.append(app)
		else:
			raise RuntimeError('Error! unknown application component name: ' + str(app.name))
			
	def remove(self, app):
		if app.name=="force":
			container = self.forces
		elif app.name=="integration":
			container = self.integrations
		elif app.name=="tinker":
			container = self.tinkers
		elif app.name=="dump":
			container = self.dumps
		else:
			raise RuntimeError('Error! unknown application component name: ' + str(app.name))

		try:
			container.remove(app)   # python list remove()
		except ValueError:
			pass
					
	def _format_hms(self, seconds):
		try:
			seconds = float(seconds)
		except Exception:
			seconds = 0.0
		if not math.isfinite(seconds) or seconds < 0.0:
			seconds = 0.0
		hours = int(seconds / 3600.0)
		seconds -= float(hours) * 3600.0
		minutes = int(seconds / 60.0)
		seconds -= float(minutes) * 60.0
		return "%02d:%02d:%02d" % (hours, minutes, int(seconds))

	def _tps_sync_if_requested(self):
		# Default ON for CuPy-native runs: make the internal application TPS a real
		# GPU elapsed-time metric instead of a kernel-launch enqueue metric.  The
		# aggregation/window rule still mirrors GALAMOST.  Disable only for low-level
		# diagnostics with PYGAMD_TPS_SYNC=0.
		flag = str(os.environ.get("PYGAMD_TPS_SYNC", "1")).strip().lower()
		if flag in ("0", "false", "no", "off"):
			return
		try:
			from pygamd.backend.cupy import import_cupy
			cp = import_cupy()
			cp.cuda.Stream.null.synchronize()
		except Exception:
			# Preserve legacy behavior if CuPy is unavailable or synchronization fails
			# during non-CuPy/static smoke tests.
			return

	def tps_compute(self, timestep):
		if timestep==self.target:
			self._tps_sync_if_requested()
			self.end_time = time.perf_counter()
			elapsed = self.end_time - self.start_time
			if elapsed <0.00000001 or elapsed>10000000.0:
				raise RuntimeError('Error! tps abnormal with time ', elapsed)			 
			tps = float(self.target-self.last) / elapsed
			all = self.end - timestep
			rseconds = float(all)/tps
			sample = {
				"tps": tps,
				"timestep": int(timestep),
				"elapsed_s": elapsed,
				"included_in_average": (not self.first_tps),
				"period_mode": self.tps_period_mode,
				"period_steps": int(self.tps_period_steps),
			}
			self.tps_samples.append(sample)
			if not self.first_tps:
				print("INFO : | TPS %5.3f"%tps, " | Time step %d"%timestep, " | Remaining time %s"%self._format_hms(rseconds))
				self.average_tps += tps
				self.count_tps += 1
				self.average_tps_value = self.average_tps/float(self.count_tps)
				self.tps_window_status = "OK_AVERAGE_TPS"
			self.last = self.target
			if self.tps_period_steps > 0:
				self.target += int(self.tps_period_steps)
			else:
				self.target += int(20.0*tps)+1
			self.start_time = self.end_time
			if self.target>self.end and self.count_tps!=0:
				print("INFO : | Average TPS : %5.3f"%(self.average_tps/float(self.count_tps)))
			self.first_tps = False

	def register(self, timestep):
		for key in self.info.compute_properties:
			self.info.compute_properties[key] = False

		for ii in self.integrations:
			ii.data.register(timestep)

		for di in self.dumps:
			di.register(timestep)

	# --- global zero-force (non-fused path): CuPy array fill, no numba --------
	def _zero_force_vp(self):
		'''Zero the global force and virial/potential buffers (non-fused, multi-
		force path) using a CuPy fill on the CuInfo device arrays.  The fused
		single-RawPair path skips this (the pair kernel overwrites force).'''
		cu = getattr(self.info, "cu", None)
		if cu is None:
			raise RuntimeError(
				"CuPy-native current implementation requires a CuInfo snapshot "
				"(pygamd.snapshot.read(..., backend='cupy')) for the global zero-force step.")
		if _sp.ENABLED: _sp.inc("force_clear_count")
		cu.d_force.fill(0.0)
		cu.d_virial_potential.fill(0.0)

	# --- fused zero-force helpers -------------------------------------
	def _fuse_target(self):
		'''Return the pair backend to fuse zero-force into, or None.

		Conservative safety conditions:
		  * exactly one force object this run, and
		  * it is a RawKernel pair force that declares can_fuse_zero_force and
		    can run in overwrite mode (set_fuse_zero_force).
		The application sees the force wrapper; the actual backend is wrapper.data.
		'''
		if len(self.forces) != 1:
			return None
		f = self.forces[0]
		target = getattr(f, "data", f)
		if bool(getattr(target, "can_fuse_zero_force", False)) and hasattr(target, "set_fuse_zero_force"):
			return target
		return None

	def _apply_fuse_zero_force(self):
		'''If safe, arm the single pair force to overwrite force this step and
		return True (caller then SKIPS the global zero_force_vp).  Otherwise put
		any RawKernel pair force back into additive mode and return False.'''
		target = self._fuse_target()
		if target is not None:
			target.set_fuse_zero_force(True)
			return True
		for fi in self.forces:
			t = getattr(fi, "data", fi)
			if hasattr(t, "set_fuse_zero_force"):
				t.set_fuse_zero_force(False)
		return False