Gamma-Ray Emission Probabilities
================================

The :math:`\gamma`-ray transmission coefficients are obtained using the strength function formalism from the expression: 

.. math::

	T^{Xl}(\epsilon_\gamma) = 2\pi f_{Xl}(\epsilon_\gamma)\epsilon_\gamma^{2l+1},

where :math:`\epsilon_\gamma` is the energy of the emitted gamma ray, :math:`Xl` is the multipolarity of the gamma ray, and :math:`f_{Xl}(\epsilon_\gamma)` is the energy-dependent gamma-ray strength function.

For :math:`E1` transitions, the `Kopecky-Uhl <https://journals.aps.org/prc/abstract/10.1103/PhysRevC.41.1941>`_ generalized Lorentzian form for the strength function is used:

.. math::

	f_{E1}(\epsilon_\gamma,T) = K_{E1}\left[ \frac{\epsilon_\gamma \Gamma_{E1}(\epsilon_\gamma)}{\left( \epsilon_\gamma^2-E_{E1}^2\right)^2 + \epsilon^2_\gamma\Gamma_{E1}(\epsilon_\gamma)^2} +\frac{0.7\Gamma_{E1}4\pi^2T^2}{E_{E1}^5} \right] \sigma_{E1}\Gamma_{E1}


where :math:`\sigma_{E1}`, :math:`\Gamma_{E1}`, and :math:`E_{E1}` are the standard giant dipole resonance (GDR) parameters. :math:`\Gamma_{E1}(\epsilon_\gamma)` is an energy-dependent damping width given by

.. math::

	\Gamma_{E1}(\epsilon_\gamma) = \Gamma\frac{\epsilon_\gamma^2+4\pi^2T^2}{E_{E1}^2},

and :math:`T` is the nuclear temperature given by

.. math::

	T=\sqrt{\frac{E^*-\epsilon_\gamma}{a(S_n)}}.

The quantity :math:`S_n` is the neutron separation energy, :math:`E^*` is the excitation energy of the nucleus, and :math:`a` is the level density parameter. The quantity :math:`K_{E1}` is obtained from normalization to experimental data on :math:`2\pi\langle \Gamma_{\gamma_0} \rangle / \langle D_0 \rangle`. 

For :math:`E2` and :math:`M1` transitions, the Brink-Axel `(Brink,1955) <https://ac.els-cdn.com/0029558287900216/1-s2.0-0029558287900216-main.pdf?_tid=8e992edb-e703-47a3-a9d3-b3182d856e85&acdnat=1520468445_d60290157ff74b45fdb22cc7519014ca>`_ `(Axel,1962) <http://journals.aps.org/pr/abstract/10.1103/PhysRev.126.671>`_ standard Lorentzian is used instead:

.. math::

	f_{Xl}(\epsilon_\gamma)=K_{Xl}\frac{\sigma_{Xl}\epsilon_\gamma\Gamma_{Xl}^2}{(\epsilon_\gamma^2-E_{Xl}^2)^2+\epsilon_\gamma^2\Gamma_{Xl}^2}.


In the current version of :program:`CGMF` (ver. |version|), only :math:`E1, E2`, and :math:`M1` transitions are allowed, and higher multipolarity transitions are neglected.


