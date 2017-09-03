classdef CoolPropWrapperTests < matlab.unittest.TestCase
% COOLPROPWRAPPERTESTS a class for testing the local system configuration
%   Please execute these tests to make sure using:  run(CoolPropWrapperTests);
%   Requires MATLAB R2014b or newer.
%   
% Copyright (C) 2017 Iliya Romm, under the MIT license.

    properties
      TestInputs2
      TestInputs6
    end
 
    methods (TestClassSetup)
      function createTestInputs(testCase)
        testCase.TestInputs2 = {'Tcrit','Water'};
        
        tmp = cell(6,1);
        tmp{1} = {'V','D'};      % Desired outputs
        tmp{2} = 'T';            % 1st input name
        tmp{3} = 200:10:6000;    % 1st input value(s)
        tmp{4} = 'P';            % 2nd input name
        tmp{5} = 5000:100:20000; % 2nd input value(s)
        tmp{6} = 'Argon';        % Fluid name          
        testCase.TestInputs6 = tmp;
      end % createTestInputs
    end % TestClassSetup methods
  
  methods (Test)
    %% Compatibility tests:
    function testPythonSetup(testCase)
      testCase.fatalAssertNotEmpty(pyversion);
      testCase.fatalAssertClass(py.importlib.import_module('numpy'),'py.module');
      % If the above fails, try reinstalling numpy.
      testCase.fatalAssertClass(py.importlib.import_module('CoolProp.CoolProp'),'py.module');
      % If the above fails, try installing the nightly release from:
      % http://www.coolprop.org/dev/coolprop/wrappers/Python/index.html#automatic-installation
    end % testPythonSetup
    
    %% Smoke tests:
    function testPythonInterfaceWithTwoInputs(testCase)
      py.CoolProp.CoolProp.PropsSI(testCase.TestInputs2{:});
    end
    
    function testPythonInterfaceWithSixInputs(testCase)
      inputs = CoolPropWrapperTests.cell2pyargs(testCase.TestInputs6);
      py.CoolProp.CoolProp.PropsSI(inputs{:});     
    end
        
    function testInterfaceEquivalence(testCase)
      SZ = [numel(testCase.TestInputs6{5}), numel(testCase.TestInputs6{3})];
      inputs = CoolPropWrapperTests.cell2pyargs(testCase.TestInputs6);      
      %% Calling PropsSI, converting python output to MATLAB and reshaping: 
      % W/o importing the CP module:
      props_ni = reshape(...
        matpy.nparray2mat(py.CoolProp.CoolProp.PropsSI(inputs{:})),...
        [SZ inputs{1}.length]);
      % W/ importing the CP module:
      CP = py.importlib.import_module('CoolProp.CoolProp');
      props_i = reshape(matpy.nparray2mat(CP.PropsSI(inputs{:})),...
        [SZ inputs{1}.length]);
      testCase.verifyEqual(props_ni,props_i);
    end % testInterfaceEquivalence
  
    %% Unit tests:
    function testHighLevelWrapperWithTwoInputs(testCase)
      tmp = testCase.TestInputs2;
      r{1} = PropsSI(tmp{:});
      r{2} = PropsSI(tmp{1},'',0,'',0,tmp{2});      
      r{3} = PropsSI({tmp{1}},tmp{2});
      testCase.verifyTrue(isequal(r{:}));
      testCase.verifyNumElements(r{1},1);
    end
    
    function testHighLevelWrapperWithSixInputs(testCase)
      tmp = testCase.TestInputs6;
      testCase.verifyNumElements(PropsSI(tmp{:}), ...
        prod( cellfun(@numel,tmp([1,3,5])) )...
      );
    end
    
    function testLowLevelWrapper(~)     
      [HEOS,CoolProp] = AbstractState('HEOS', 'Water');
      % Set a very low density state point
      HEOS.update(CoolProp.DmolarT_INPUTS, 1e-6, 300);      
      % Specify the phase:
      HEOS.specify_phase(CoolProp.iphase_gas);
      % Specify delimited fluid names:
      HEOS = AbstractState('HEOS', 'Methane&Ethane');
      % Set the mole fractions of the mixture:
      HEOS.set_mole_fractions([0.2,0.8]);
      % Do the dewpoint calculation:
      HEOS.update(CoolProp.PQ_INPUTS, 101325, 1);
      % Liquid phase molar density:
      HEOS.saturated_liquid_keyed_output(CoolProp.iDmolar);
      % Vapor phase molar density:
      HEOS.saturated_vapor_keyed_output(CoolProp.iDmolar);
      % Liquid phase mole fractions:
      HEOS.mole_fractions_liquid();
      % Vapor phase mole fractions:
      HEOS.mole_fractions_vapor();
    end
    
  end % Test methods
  
  methods (Access = private, Static = true)
    function pyInputs = cell2pyargs(input)
      assert(numel(input) == 6);
      [XX,YY] = meshgrid(input{3}, input{5}); 
      pyInputs = {py.list(input{1}),...
        input{2}, py.numpy.ravel(XX(:).'), ...
        input{4}, py.numpy.ravel(YY(:).'), ...
        input{6}};
    end
  end % private static methods
end