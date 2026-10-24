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

# -----------------------------------------------------------------------------
# Central opt-in profiler (profiler only; no physics / API change).
#
# Disabled by default.  Enable with the environment variable:
#       PYGAMD_PROFILE=1
# Optional:
#       PYGAMD_PROFILE_PERIOD=1000        (steps between stdout summaries)
#       PYGAMD_PROFILE_JSON=path.json     (final JSON dump)
#       PYGAMD_PROFILE_TEXT=path.txt      (final text dump)
#
# When disabled, every entry point is a single attribute read + early return, so
# the RawKernel-native main path runs essentially unchanged.
#
# IMPORTANT (reported to the user): when enabled, GPU stages are measured as
# *synchronized wall time* -- begin()/end() call a device synchronize so the
# elapsed time reflects the GPU work of that stage.  This adds synchronization
# overhead and inflates the absolute step time, so profiler-enabled numbers are
# meant for locating proportions, NOT as production TPS.  cupy is imported
# lazily; `import pygamd` never depends on cupy or on this profiler being on.
# -----------------------------------------------------------------------------

import os
import time
import json


def _truthy(v):
	return v not in ("", "0", "false", "False", "no", "off", None)


class Profiler:
	def __init__(self, enabled=False, period=1000, json_path=None, text_path=None):
		self.enabled = bool(enabled)
		self.period = int(period) if period else 1000
		self.json_path = json_path
		self.text_path = text_path

		# name -> [calls, total_ms];  name -> int
		self._times = {}
		self._counts = {}
		self._starts = {}
		self._step = 0
		self._cp = None          # lazy cupy handle (None=unknown, False=absent)
		self._sync_calls = 0

	# --- construction from environment ---------------------------------------
	@classmethod
	def from_env(cls):
		enabled = _truthy(os.environ.get("PYGAMD_PROFILE", ""))
		period = os.environ.get("PYGAMD_PROFILE_PERIOD", "1000")
		try:
			period = int(period)
		except Exception:
			period = 1000
		json_path = os.environ.get("PYGAMD_PROFILE_JSON", None)
		text_path = os.environ.get("PYGAMD_PROFILE_TEXT", None)
		return cls(enabled=enabled, period=period, json_path=json_path, text_path=text_path)

	# --- device synchronization (lazy cupy only; no numba) -------------------
	def _get_cp(self):
		if self._cp is False:
			return None
		if self._cp is None:
			try:
				import cupy as cp
				self._cp = cp
			except Exception:
				self._cp = False
				return None
		return self._cp

	def _sync(self):
		self._sync_calls += 1
		cp = self._get_cp()
		if cp is not None:
			try:
				cp.cuda.runtime.deviceSynchronize()
			except Exception:
				pass

	# --- timing API (synchronized wall time) ---------------------------------
	def begin(self, name):
		if not self.enabled:
			return
		self._sync()
		self._starts[name] = time.perf_counter()

	def end(self, name):
		if not self.enabled:
			return
		self._sync()
		t0 = self._starts.get(name)
		if t0 is None:
			return
		dt_ms = (time.perf_counter() - t0) * 1000.0
		e = self._times.get(name)
		if e is None:
			self._times[name] = [1, dt_ms]
		else:
			e[0] += 1
			e[1] += dt_ms

	# allow recording an externally measured time (e.g. CUDA-event ms)
	def add_time(self, name, ms, calls=1):
		if not self.enabled:
			return
		e = self._times.get(name)
		if e is None:
			self._times[name] = [calls, float(ms)]
		else:
			e[0] += calls
			e[1] += float(ms)

	def count(self, name, inc=1):
		if not self.enabled:
			return
		self._counts[name] = self._counts.get(name, 0) + inc

	# --- per-step hook --------------------------------------------------------
	def tick_step(self):
		if not self.enabled:
			return
		self._step += 1
		if self.period > 0 and self._step % self.period == 0:
			self.report()

	# --- reporting ------------------------------------------------------------
	def _steps(self):
		return self._step if self._step > 0 else 1

	def report(self, step=None):
		if not self.enabled:
			return
		steps = self._steps()
		lines = []
		lines.append("[pygamd profile] step=%d steps=%d (synchronized wall time)" % (self._step, steps))
		# sort items by total_ms descending
		for name, (calls, total_ms) in sorted(self._times.items(), key=lambda kv: -kv[1][1]):
			lines.append("  %-32s %8.4f ms/step   (calls=%d, total=%.2f ms)"
			             % (name, total_ms / steps, calls, total_ms))
		if self._counts:
			lines.append("  -- counters --")
			for name, c in sorted(self._counts.items()):
				lines.append("  %-32s %d" % (name, c))
		print("\n".join(lines))

	def _as_dict(self):
		steps = self._steps()
		items = {}
		for name, (calls, total_ms) in self._times.items():
			items[name] = {
				"calls": calls,
				"total_ms": round(total_ms, 4),
				"avg_ms": round(total_ms / calls, 6) if calls else 0.0,
				"ms_per_step": round(total_ms / steps, 6),
			}
		return {"steps": self._step, "sync_calls": self._sync_calls,
		        "items": items, "counters": dict(self._counts)}

	def write_json(self, path=None):
		if not self.enabled:
			return
		path = path or self.json_path
		if not path:
			return
		try:
			with open(path, "w") as f:
				json.dump(self._as_dict(), f, indent=2)
		except Exception as e:
			print("info : [pygamd profile] failed to write JSON: %s" % repr(e))

	def write_text(self, path=None):
		if not self.enabled:
			return
		path = path or self.text_path
		if not path:
			return
		steps = self._steps()
		try:
			with open(path, "w") as f:
				f.write("pygamd profile (synchronized wall time)\n")
				f.write("steps=%d  sync_calls=%d\n\n" % (self._step, self._sync_calls))
				f.write("--- timings (sorted by total_ms) ---\n")
				for name, (calls, total_ms) in sorted(self._times.items(), key=lambda kv: -kv[1][1]):
					f.write("%-34s total=%10.3f ms   ms/step=%8.5f   avg=%9.6f ms   calls=%d\n"
					        % (name, total_ms, total_ms / steps, (total_ms / calls if calls else 0.0), calls))
				f.write("\n--- timings (sorted by ms/step) ---\n")
				for name, (calls, total_ms) in sorted(self._times.items(), key=lambda kv: -kv[1][1] / steps):
					f.write("%-34s ms/step=%8.5f   total=%10.3f ms\n" % (name, total_ms / steps, total_ms))
				f.write("\n--- counters ---\n")
				for name, c in sorted(self._counts.items()):
					f.write("%-34s %d\n" % (name, c))
		except Exception as e:
			print("info : [pygamd profile] failed to write TEXT: %s" % repr(e))

	def finalize(self):
		if not self.enabled:
			return
		print("info : [pygamd profile] === final summary ===")
		self.report()
		self.write_json()
		self.write_text()


# module-level singleton, configured from the environment at import time
prof = Profiler.from_env()
