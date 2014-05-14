Forming the Phase Boundary
==========================

Overview
--------
The analysis in this section follows the methodologies proposed in the GERG 2004 monograph published in 2007

System of Equations
-------------------

Our residual vector :math:`\mathbf{F}` is equal to 

.. math::

    F_i = \ln\phi(T,p,\mathbf{y})-\ln \phi(T,p,\mathbf{x})+\ln K_i=0,  i=1,2,3... N

.. math::

    F_{N+1} = \sum_{i=1}^{N}(y_i-x_i)=0
    
.. math::

    x_i = \frac{z_i}{1-\beta+\beta K_i}
    
and 

.. math::

    y_i = \frac{K_iz_i}{1-\beta+\beta K_i}
    
.. math::

    F_{N+2} = X_s - S = 0
    
and the system to be solved is equal to

.. math::

    \mathbf{J}\mathbf{\Delta X}= -\mathbf{F}
    

Building the Jacobian matrix
----------------------------
This is the trickiest part of this method.  There are a lot of derivatives to implement, and we want to implement all of them analytically.

.. math::

    \frac{\partial F_i}{\partial \ln T} = T\left[ \left(\frac{\partial \ln \phi_i}{\partial T}\right)''_{p,\mathbf{n}} -\left(\frac{\partial \ln \phi_i}{\partial T}\right)'_{p,\mathbf{n}}\right]
    
.. math::

    \frac{\partial F_i}{\partial \ln p} = p\left[ \left(\frac{\partial \ln \phi_i}{\partial p}\right)''_{T,\mathbf{n}} -\left(\frac{\partial \ln \phi_i}{\partial p}\right)'_{T,\mathbf{n}}\right]
    
.. math::

    \frac{\partial F_i}{\partial \ln K_j} = \frac{K_jz_j}{(1-\beta+\beta K_j)^2}[(1-\beta)\phi_{ij}''+\beta\phi_{ij}']+\zeta

where :math:`\zeta = 0` for i:math:`\neq`j , and  :math:`\zeta = 0` for i=j.  Also

.. math::

    \phi_{ij} = n\left( \frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p}

For the :math:`F_{N+1}` term,

.. math::

    \frac{\partial F_{N+1}}{\partial \ln K_j}=\frac{K_jz_j}{(1-\beta+\beta K_j)^2}    

and all other partials of :math:`F_{N+1}` in the Jacobian are zero.  For the specified term

.. math::

    \frac{\partial F_{N+2}}{X_s}=1
    
and all other partials of :math:`F_{N+2}` in the Jacobian are zero.
    
..

    Onwards...

    Gerg 2004 Monograph, Eqn 7.27:

    .. math::

        \ln \phi_i  = \left( \frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}-\ln Z
        
    and (Kunz, 2012, Table B4)

    .. math::

        \left( \frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j} = \alpha^r + n\left( \frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}
        
    so

    .. math::

        \ln \phi_i  = \alpha^r + n\left( \frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}-\ln Z
    
    and its derivative w.r.t T can be obtained analytically.  What about pressure?


The fugacity coefficient can be obtained from (Kunz, 2012, equation 29)

From GERG Monograph p. 60: 

    Since the two phases of a non-critical mixture are characterised by different compositions resulting in different values for the reducing functions and the corresponding reduced variables, a simple integral criterion which connects all phase equilibrium properties in a single relation such as Eq. (4.11) does not exist for mixtures

Pandoc
------

pandoc --mathjax -s -f rst -t html5 -o phase_boundary.html phase_boundary.rst