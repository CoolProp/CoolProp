function props = PropsSI(varargin)
% PROPSSI is a MATLAB wrapper for the python binding of the CoolProp 
%   high-level interface - PropsSI.
%
%   The purpose of this tool is to allow convenient and efficient (vectorized) 
%   evaluation of PropsSI from MATLAB by leveraging the well-maintained and 
%   documented python binding of CoolProp.PropsSI()
%   
% INPUTS:
% Number of supported inputs: 0 (demo case), 2 (trivial properties), 6 (all the rest). 
%
% varargin{1} - *char array* or a *cell array of char arrays* representing the 
%   desired outputs, according to:
%   http://www.coolprop.org/coolprop/HighLevelAPI.html#table-of-string-inputs-to-propssi-function
% varargin{2} - A *char array* representing the meaning of the next input vector.
% varargin{3} - A vector of doubles.
% varargin{4} - A *char array* representing the meaning of the next input vector.
% varargin{5} - A vector of doubles.
% varargin{6} - a *char array* representing the fluid name, according to:
%   http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
% Specifying a mixture is possible, but requires additional steps as shown in:
%   http://www.coolprop.org/fluid_properties/Mixtures.html
% 
% Copyright (C) 2017 Iliya Romm, under the MIT license.

%% Handling inputs:
switch nargin
  case 0 % Demonstration case
  % Below is a demonstration of calling the 6-argument variant of PropsSI().
  % Besides being an example of usage, it also verifies that the MATLAB->python 
  %  interface was set up correctly.
    demoCase = cell(6,1);
    demoCase{1} = {'V','D'};      % Desired outputs
    demoCase{2} = 'T';            % 1st input name
    demoCase{3} = 200:10:6000;    % 1st input value(s)
    demoCase{4} = 'P';            % 2nd input name
    demoCase{5} = 5000:100:20000; % 2nd input value(s)
    demoCase{6} = 'Argon';        % Fluid name
    props = PropsSI(demoCase{:});
    return
  case 2 % Convert the 2-parameter call into a 6-parameter one:
    toSixAdapter = cell(6,1);
    toSixAdapter{1} = varargin{1};
    toSixAdapter{2} = '';
    toSixAdapter{3} = 0;
    toSixAdapter{4} = '';
    toSixAdapter{5} = 0;
    toSixAdapter{6} = varargin{2};    
    props = PropsSI(toSixAdapter{:});
    return    
  case 6 % If correct # of inputs, perform input validation:
    validateattributes(varargin{1}, {'char','cell'},{'nonempty'});
    if iscell(varargin{1})
      multiOutputs = true;
      assert(all(cellfun(@ischar,varargin{1})));
    else
      multiOutputs = false;
    end
    assert(ischar(varargin{2}));
    validateattributes(varargin{3}, {'numeric'},  {'vector'});
    assert(ischar(varargin{4}));
    validateattributes(varargin{5}, {'numeric'},  {'vector'});
    assert(ischar(varargin{6}));
  otherwise
    throw(MException('PropsSI:badNumInputs',...
                     'Unsupported number of inputs!'));
end
%% Converting varargin to python-friendly inputs:
[XX,YY] = meshgrid(varargin{3}, varargin{5});
if multiOutputs
  inputs = {py.list(varargin{1}), varargin{2}, py.numpy.ravel(XX(:).'), ...
                                  varargin{4}, py.numpy.ravel(YY(:).'), ...
                                  varargin{6}};
else
  inputs = {varargin{1}, varargin{2}, py.numpy.ravel(XX(:).'), ...
                         varargin{4}, py.numpy.ravel(YY(:).'), ...
                         varargin{6}};
end
%% Calling PropsSI, converting python output to MATLAB and reshaping: 
if multiOutputs
  props = reshape(matpy.nparray2mat(py.CoolProp.CoolProp.PropsSI(inputs{:})),...
                  [size(XX) inputs{1}.length]);
else
  props = reshape(matpy.nparray2mat(py.CoolProp.CoolProp.PropsSI(inputs{:})),...
                  [size(XX) 1]);      
end