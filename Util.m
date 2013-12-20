classdef Util<handle
    
    methods(Static)
        function Aout = displace(A, xt, yt)
            % Displaces input array A by 'xt' horizontally and 'yt' vertically, by
            % adding zeros and discarding elements falling out of the original array
            % size.
            %
            % SYNTAX:
            % Aout = Displace(Ain, xt, yt)
            %
            % DESCRIPTION:
            % If inputs 'xt' and 'yt' are scalar (necessarily integers) and xt, yt > 0,
            % the function output will be input matrix displace to the right and down,
            % respectively.
            %
            % If 'Ain' has dimension higher than 2 and 'xt' and 'yt' are scalars, the
            % same displacement operation is performed on each matrix (corresponding to
            % the first two dimensions of 'Ain') from the set of all matrices (each
            % associated with the indexes from remaining dimensions). That is, this is
            % equivalent to performing Displace(Ain(:,:,i), xt, yt) for all i.
            %
            % If 'Ain' has dimension higher than 2 and 'xt' and/or 'yt' are arrays (be
            % it matrices or vectors), each matrix that forms 'Ain' is displaced an
            % amount associated with the corresponding elements of 'xt' and 'yt'. That
            % is, this is equivalent to performing Displace(Ain(:,:,i), xt(i), yt(i))
            % for all i. A(1,1,:), xt(:) and yt(:) should have the same length or one
            % of the latter two can still be a scalar. All elements of 'xt' and 'yt'
            % should be integers.
            
            % Input argument error checking.
            if nargin ~= 3
                error('displace:argChk', 'Wrong number of input arguments. It should be 3.')
            end
            if ~all(round(xt) == xt) || ~all(round(yt) == yt)
                error('displace:intChk', 'Second and third arguments must be integers.')
            end
            
            % Vector with sizes of A
            m = zeros(1,ndims(A));
            for i = 1 : ndims(A)
                m(i) = size(A,i);
            end
            
            % Number of matrices forming A
            len = length(A(1,1,:));
            
            % Extend input displacements if scalars
            if isscalar(xt)
                xt = ones(len,1)*xt;
            end
            if isscalar(yt)
                yt = ones(len,1)*yt;
            end
            
            % Check length match between input matrix and displacement arrays
            if (length(xt(:)) ~= len) || (length(yt(:)) ~= len)
                error('displace:argSizeChk', ...
                    'Size mismatch between input matrix and displacement vector.')
            end
            
            % Vectors with sizes of displacement matrices in X and Y directions
            mx = [repmat(m(1),len,1), abs(xt(:))];
            my = [abs(yt(:)), repmat(m(2),len,1)];
            
            
            Aout = A;
            
            % Horizontal displacement
            for i = 1 : len
                if xt(i) ~= 0
                    if xt(i) > 0
                        Atemp = [zeros(mx(i,:)), Aout(:,:,i)];
                        Aout(:,:,i) = Atemp(:, 1: end - xt(i));
                    elseif xt(i) < 0
                        Atemp = [ Aout(:,:,i), zeros(mx(i,:))];
                        Aout(:,:,i) = Atemp(:, 1 + abs(xt(i)) : end);
                    end
                end
            end
            Aout = reshape(Aout, m);
            
            % Vertical displacement
            for i = 1 : len
                if yt(i) ~= 0
                    if yt(i) > 0
                        Atemp = [zeros(my(i,:)); Aout(:,:,i)];
                        Aout(:,:,i) = Atemp(1: end - yt(i), :);
                    elseif yt(i) < 0
                        Atemp = [ Aout(:,:,i); zeros(my(i,:))];
                        Aout(:,:,i) = Atemp(1 + abs(yt(i)): end, :);
                    end
                end
            end
            Aout = reshape(Aout, m);
        end
        function Aout = crop(A, width, height)
            if nargin ~= 3
                error('crop:argChk', 'Wrong number of input arguments. It should be 3.')
            end
            if ~all(round(width) == width) || ~all(round(height) == height)
                error('crop:intChk', 'Second and third arguments must be integers.')
            end
            
            % Vector with sizes of A
            m = zeros(1,ndims(A));
            for i = 1 : ndims(A)
                m(i) = size(A,i);
            end
            
            if (width == m(2) || height == m(1))
                Aout = A;
                return;
            end
            
            centerIndex = [floor(m(1)/2 + 1/2), floor(m(2)/2+1/2)];
            Aout = A(centerIndex(1) - floor(height/2) : ...
                centerIndex(1) - floor(height/2) + height - 1, ...
                centerIndex(2) - floor(width/2) : ...
                centerIndex(2) - floor(width/2) + width - 1, :);
        end
        function str = transformInputParametersIntoStructure(params)
            str = struct;
            for i = 1 : floor(numel(params)/2)
                 str.(params{2*i-1}) = params{2*i};
            end
        end
        function normArray = normalize(A)
        % normalize Normalize matrices that compose n-dimensional input array
            
            normArray = A;
            numberOfMatrices = length(A(1,1,:));
            for i = 1 : numberOfMatrices
                normArray(:,:,i) = A(:,:,i)/sum(sum(A(:,:,i)));
            end          
        end
    end
end