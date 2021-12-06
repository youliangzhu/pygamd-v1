.. _page-units:

Units
=====

PYGAMD stores and computes all values in reduced units. The quantities in real units can be converted into 
the ones in reduced units by defining a set of fundamental units by user himself.

Fundamental Units
-----------------

The three fundamental units are:

- distance - :math:`\mathcal{\sigma}`
- energy - :math:`\mathcal{\varepsilon}`
- mass - :math:`\mathcal{m}`


Temperature units (thermal energy)
----------------------------------

PYGAMD accepts all temperature inputs and provides all temperature output values in units of energy:
:math:`k_{B} T`, where :math:`k_{B}` is Boltzmann's constant. In reduced units, one usually reports the value
:math:`T^* = k_{B}T/\mathcal{\varepsilon}`.

.. _charge-units:

Charge units
------------

The charge used in PYGAMD is also reduced. The units of charge are: 
:math:`(4 \pi \epsilon_0 \epsilon_r \mathcal{\sigma} \mathcal{\varepsilon})^{1/2}`, where
:math:`\epsilon_0` is vacuum permittivity and :math:`\epsilon_r` is relative permittivity.

With :math:`f= 1/4\pi \epsilon_0=138.935\text{ }kJ\text{ }mol^{-1}\text{ }nm\text{ }e^{-2}`,
the units of charge are: 
:math:`(\epsilon_r \mathcal{\sigma} \mathcal{\varepsilon}/f)^{1/2}`.
Divide a given charge by this quantity to convert it into an input value for PYGAMD.

Common derived units
--------------------

Here are some commonly used derived units:

- time - :math:`\tau = \sqrt{\mathcal{m} \mathcal{\sigma}^2/\mathcal{\varepsilon}}`
- volume - :math:`\mathcal{\sigma}^3`
- velocity - :math:`\mathcal{\sigma}/\tau`
- momentum - :math:`\mathcal{m}\mathcal{\sigma}/\tau`
- acceleration - :math:`\mathcal{\sigma}/\tau^2`
- force - :math:`\mathcal{\varepsilon}/\mathcal{\sigma}`
- pressure - :math:`\mathcal{\varepsilon}/\mathcal{\sigma}^3`

Example physical units
----------------------

There are many possible choices of physical units that one can assign. One common choice is:

- distance - :math:`\mathcal{\sigma} = \mathrm{nm}`
- energy - :math:`\mathcal{\varepsilon} = \mathrm{kJ/mol}`
- mass - :math:`\mathcal{m} = \mathrm{amu}`

Derived units / values in this system:

- time - picoseconds
- velocity - nm/picosecond
- pressure - 16.3882449645417 atm
- force - 1.66053892103218 pN
- :math:`k_{B}` = 0.00831445986144858 kJ/mol/Kelvin
