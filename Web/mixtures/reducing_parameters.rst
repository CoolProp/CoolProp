
Mixtures Information
====================


Reducing Parameters
-------------------

From Lemmon, JPCRD, 2000 for the properties of Dry Air, and also from Lemmon, JPCRD, 2004 for the properties of R404A, R410A, etc.

.. math::

    \rho_r(\bar x) = \left[ \sum_{i=1}^m\frac{x_i}{\rho_{c_i}}+\sum_{i=1}^{m-1}\sum_{j=i+1}^{m}x_ix_j\zeta_{ij}\right]^{-1}

.. math::

    T_r(\bar x) = \sum_{i=1}^mx_iT_{c_i}+\sum_{i=1}^{m-1}\sum_{j=i+1}^mx_ix_j\xi_{ij}

From the GERG 2008 formulation (Kunz and Wagner, JCED, 2012)

.. math::

    T_r(\bar x) = \sum_{i=1}^{N}x_i^2T_{c,i} + \sum_{i=1}^{N-1}\sum_{j=i+1}^{N}2x_ix_j\beta_{T,ij}\gamma_{T,ij}\frac{x_i+x_j}{\beta_{T,ij}^2x_i+x_j}(T_{c,i}T_{c,j})^{0.5}
    
.. math::

    \frac{1}{\rho_r(\bar x)}=v_r(\bar x) = \sum_{i=1}^{N}x_i^2\frac{1}{\rho_{c,i}} + \sum_{i=1}^{N-1}\sum_{j=i+1}^N2x_ix_j\beta_{v,ij}\gamma_{v,ij}\frac{x_i+x_j}{\beta^2_{v,ij}x_i+x_j}\frac{1}{8}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}
    
Excess Helmholtz Energy Terms
-----------------------------
From Lemmon, JPCRD, 2004 for the properties of R404A, R410A, etc.

.. math::

    \alpha^E(\delta,\tau,\mathbf{x}) = \sum_{i=1}^{m-1} \sum_{j=i+1}^{m} \left [ x_ix_jF_{ij} \sum_{k}N_k\delta_{d_k}\tau^{t_k}\exp(-\delta^{l_k})\right]
    
where the terms :math:`N_k,d_k,t_k,l_k` correspond to the pair given by the indices :math:`i,j`

From Lemmon, JPCRD, 2000 for the properties of Dry Air

.. math::

    \alpha^E(\delta,\tau,\mathbf{x}) = \left \lbrace \sum_{i=1}^{2} \sum_{j=i+1}^{3} x_ix_jF_{ij}\right\rbrace \left[-0.00195245\delta^2\tau^{-1.4}+0.00871334\delta^2\tau^{1.5} \right]


From Kunz and Wagner, JCED, 2012 for GERG 2008

.. math::

    \alpha^E(\delta,\tau,\mathbf{x}) = \sum_{i=1}^{N-1} \sum_{j=i+1}^{N} x_ix_jF_{ij}\alpha_{ij}^r(\delta,\tau)
    
where

.. math::

    \alpha_{ij}^r(\delta,\tau) = \sum_{k=1}^{K_{pol,ij}}\eta_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}}+\sum_{k=K_{pol,ij}+1}^{K_{pol,ij}+K_{Exp,ij}}\eta_{ij,k}\delta^{d_{ij,k}}\tau^{t_{ij,k}}\exp[-\eta_{ij,k}(\delta-\varepsilon_{ij,k})^2-\beta_{ij,k}(\delta-\gamma_{ij,k})]
    
and is for the particular binary pair given by the indices :math:`i,j`.  This term is similar in form to other Helmholtz energy terms for pure fluids though the derivatives are slightly special.

Appendix
--------
To convert from the form from Lemmon for HFC and Air to that of GERG 2008, the following steps are required:

.. math::

    x_0T_{c0}+(1-x_0)T_{c1}+x_0(1-x_0)\xi_{01} = x_0^2T_{c0}+(1-x_0)^2T_{c1} + 2x_0(1-x_0)\beta\gamma_T\frac{x_0+(1-x_0)}{\beta x_0 + (1-x_0)}\sqrt{T_{c0}T_{c1}}
    
set :math:`\beta=1`, solve for :math:`\gamma`.  Equate the terms

.. math::

    x_0T_{c0}+(1-x_0)T_{c1}+x_0(1-x_0)\xi_{01} = x_0^2T_{c0}+(1-x_0)^2T_{c1} + 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}
    
Move to LHS

.. math::

    [x_0-x_0^2]T_{c0}+[(1-x_0)-(1-x_0)^2]T_{c1}+x_0(1-x_0)\xi_{01} = 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}

Factor

.. math::

    x_0(1-x_0)T_{c0}+(1-x_0)[1-(1-x_0)]T_{c1}+x_0(1-x_0)\xi_{01} = 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}
    
Expand

.. math::

    x_0(1-x_0)T_{c0}+x_0(1-x_0)T_{c1}+x_0(1-x_0)\xi_{01} = 2x_0(1-x_0)\gamma_T\sqrt{T_{c0}T_{c1}}
    
Cancel factors of :math:`x_0(1-x_0)`

.. math::

    T_{c0}+T_{c1}+\xi_{01} = 2\gamma_T\sqrt{T_{c0}T_{c1}}
    
Answer:

.. math::

    \boxed{\gamma_T = \dfrac{T_{c0}+T_{c1}+\xi_{01}}{2\sqrt{T_{c0}T_{c1}}}}
    
Same idea for the volume

.. math::

    \boxed{\gamma_v = \dfrac{v_{c0}+v_{c1}+\zeta_{01}}{\frac{1}{4}\left(\frac{1}{\rho_{c,i}^{1/3}}+\frac{1}{\rho_{c,j}^{1/3}}\right)^{3}}}