classdef matpy
  %MATPY A class for converting between MATLAB and numpy (python) ndarrays
  %
  % Based on https://www.mathworks.com/matlabcentral/answers/157347
  % by David Laurenson (2016), Christoph Wiedemann (2017), and Iliya Romm (2017).
  
  methods (Access = public, Static = true)
    
    function result = mat2nparray( mtlbArr )
      %mat2nparray Convert a MATLAB array into a numpy nparray
      %   Convert an n-dimensional MATLAB array into an equivalent nparray
      data_size = size(mtlbArr);
      switch numel(data_size)
        case {0,1}
          % scalars and 1-D vectors are trivial
          result = py.numpy.array(mtlbArr);
        case 2
          % A transpose operation is required either in MATLAB, or in Python due
          % to the difference between row major and column major ordering
          transposed = mtlbArr.';
          % Pass the array to Python as a vector, and then reshape to the correct
          % size
         result = py.numpy.reshape(transposed(:).', int32(data_size));
        otherwise
          % For an n-dimensional array, transpose the first two dimensions to
          % sort the storage ordering issue
          transposed = permute(mtlbArr, numel(data_size):-1:1);
          % Pass it to python, and then reshape to the python style of matrix
          % sizing
          result = py.numpy.reshape(transposed(:).', int32(fliplr(size(transposed))));
      end
    end
    
    function result = nparray2mat( npArr )
      %nparray2mat Convert an nparray from numpy to a MATLAB array
      %   Convert an n-dimensional nparray into an equivalent MATLAB array
      data_size = cellfun(@int64, cell(npArr.shape));
      if ~all(data_size) % Array with at least one zero dimension
        result = zeros(data_size);
        return
      end
      switch numel(data_size)
        case {0,1}
          % Convert to a double scalar or vector (preserving shape)
          result = double(py.array.array('d', py.numpy.nditer(npArr)));
        case 2
          % order='F' is used to get data in column-major order (as in Fortran
          % 'F' and MATLAB)
          result = reshape(double(py.array.array('d', ...
            py.numpy.nditer(npArr, pyargs('order', 'F')))), ...
            data_size);
        otherwise
          % For multidimensional arrays more manipulation is required
          % First recover in python order (C contiguous order)
          result = double(py.array.array('d', ...
            py.numpy.nditer(npArr, pyargs('order', 'C'))));
          % Switch the order of the dimensions (as Python views this in the
          % opposite order to MATLAB) and reshape to the corresponding C-like
          % array
          result = reshape(result, fliplr(data_size));
          % Now transpose rows and columns of the 2D sub-arrays to arrive at the
          % correct MATLAB structuring
          result = permute(result, numel(data_size):-1:1);
      end
    end
  end
end
