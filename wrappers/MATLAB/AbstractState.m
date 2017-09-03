function [abState, CoolProp] = AbstractState(varargin)
% AbstractState is a MATLAB wrapper for the python binding of the CoolProp 
%   low-level interface - AbstractState.
%
% It is highly advised to assign the 2nd output of this function into a variable 
% called 'CoolProp' in the calling workspace (for compatibility with the python 
% code examples appearing in the documentation).
%
% IMPORTANT: 
% (1) No input checking is performed! It is up to the user to ensure valid inputs.
% (2) When using the state object for computations, outputs should be assumed to
%     consist of python objects (arrays, lists, etc.) and converted to MATLAB 
%     variables accordingly. The provided matpy utility may be useful for that.
% 
% INPUTS:
%   The 1st input contains the backend name.
%   The 2nd input contains a fluid name, or a set of fluid names delimited by '&'.
%
%   The following input structures are supported:
%   AbstractState({'HEOS', 'Water'});
%   AbstractState( 'HEOS', 'Water' );
%   AbstractState( 'BICUBIC&HEOS', 'Methane&Ethane' );
%
%   In MATLAB R2016b and newer, string objects can be used in place of char arrays: 
%      AbstractState( "HEOS", "Water" ); % ( >=R2016b syntax )
% 
% For more information on CoolProp's AbstractState class see:
% 1. http://coolprop.sourceforge.net/coolprop/LowLevelAPI.html
% 2. \CoolProp\include\AbstractState.h
%
% Copyright (C) 2017 Iliya Romm, under the MIT license.

CoolProp = py.importlib.import_module('CoolProp.CoolProp');
if iscell(varargin{1})
  abState = CoolProp.AbstractState(varargin{1}{:});
else
  abState = CoolProp.AbstractState(varargin{:});
end