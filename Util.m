classdef Util
    
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
            
            % Return if displacement is set to 0
            if (all(~xt) && all(~yt))
                Aout = A;
                return;
            end

            % Vector with sizes of A
            m = zeros(1,ndims(A));
            m = size(A);
            
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
            
            Aout = zeros(size(A));
            
            % Horizontal displacement
            for i = 1 : len
                if xt(i) > 0
                    Aout(:,abs(xt(i))+1:end,i) = A(:, 1: end - abs(xt(i)), i);
                elseif xt(i) < 0
                    Aout(:, 1: end -abs(xt(i)), i) = A(:, 1 + abs(xt(i)) : end, i);
                end
            end

            if size(Aout)~=m
                Aout = reshape(Aout, m);
            end
            
            % Vertical displacement
            for i = 1 : len
                if yt(i) > 0
                    Aout(abs(yt(i))+1:end,:,i) = A(1: end - abs(yt(i)), :,i);
                elseif yt(i) < 0
                    Aout(1: end-abs(yt(i)),:,i) = A(1 + abs(yt(i)): end, :, i);
                end
            end

            if size(Aout)~=m
                Aout = reshape(Aout, m);
            end

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
            
            if (width == m(2) && height == m(1))
                Aout = A;
                return;
            end
            
            centerIndex = [ceil(m(1)/2 + 1/2), ceil(m(2)/2+1/2)];
            Aout = A(centerIndex(1) - floor(height/2) : ...
                centerIndex(1) - floor(height/2) + height - 1, ...
                centerIndex(2) - floor(width/2) : ...
                centerIndex(2) - floor(width/2) + width - 1, :);
        end
        
        function B = rot90All( A, n )
            %rot90All Rotates multidimensional array.
            %   SYNTAX
            %   B = rot90All(A,n)
            %
            %   DESCRIPTION
            %   Applicable when A has dimensionality higher than 2, otherwise
            %   reproduces behavior of function rot90.
            %
            %   When A has dimensions higher than 2, rotates matrices associated with
            %   each set of indexes that correspond to dimensions 3 and higher. For
            %   example, suppose
            %   A = rand(3,2,2,3),
            %   then rot90All will treat each matrix A(:,:,1,1), A(:,:,1,2), ...
            %   A(:,:,end,end) separately, and apply rot90 to them.
            %
            %   Output will have same size as input.
            
            dim = ndims(A);
            
            if nargin == 1
                n = 1;
            end
            
            % Vector with sizes of A
            m = zeros(1,dim);
            for i = 1 : dim
                m(i) = size(A,i);
            end
            
            len = length(A(1,1,:));
            B = zeros(m(1),m(2),len);
            
            for i = 1 : len
                B(:,:,i) = rot90(A(:,:,i),n);
            end
            
            B = reshape(B, m);
            
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

        function [optRow, optCol] = findOptimumSubplotGrid(num)
            optRow = 1;
            optCol = 1;
            while true
                if (optRow * optCol >= num)
                    return
                elseif (optCol == optRow)
                    optCol = optCol + 1;
                else
                    optRow = optRow + 1;
                end
            end
        end

        function pvtProp = getSetPrivateProperties(obj)
            pvtProp = {};
            allProp = properties(obj);
            for p = 1 : length(allProp)
                metaProp = findprop(obj, allProp{p});
                if strcmpi(metaProp.SetAccess, 'private')
                    pvtProp = [pvtProp; allProp{p}];
                end
            end
        end
        
        function pblProp = getSetPublicProperties(obj)
            pblProp = {};
            allProp = properties(obj);
            for p = 1 : length(allProp)
                metaProp = findprop(obj, allProp{p});
                if strcmpi(metaProp.SetAccess, 'public')
                    pblProp = [pblProp; allProp{p}];
                end
            end
        end
        
        function ip = modeInnerProduct(A, B)
            [a1, a2, a3] = size(A);
            [b1, b2, b3] = size(B);

            if (a1 ~= b1 || a2 ~= b2)
                error('util:modeInnerProduct', 'Input arguments should be composed of matrices that are the same size.');
            end
            ip = zeros(a3,b3);

            for ib = 1 : b3
                for ia = 1 : a3
                    ip(ia,ib) = sum(sum(conj(A(:,:,ia)) .* B(:,:,ib), 2), 1);
                    ip(ia,ib) = ip(ia,ib)/sqrt(sum(sum(abs(A(:,:,ia)).^2 ,2),1) * sum(sum(abs(B(:,:,ib)).^2 ,2),1));

                end
            end
        end

        function prt = getParityComponents(modes)
            evenComponent = modes+Util.rot90All(modes,2);
            oddComponent = modes-Util.rot90All(modes,2);

            even = sum(sum(abs(evenComponent).^2,2),1);

            odd = sum(sum(abs(oddComponent).^2,2),1);
            prt(:,1) = even(:)./(even(:)+odd(:));
            prt(:,2) = odd(:)./(even(:)+odd(:));
        end
        
        function [nx, ny] = getHermiteGaussOrders(num)
            nx = floor(num/10);
            ny = floor(num - 10*nx);
        end

        function [nxMax, nyMax] = getMaximumHermiteGaussOrders(m)
            nOrd = numel(m);
            nxMax = 0;
            nyMax = 0;

            for i = 1 : nOrd
                [nx,ny] = Util.getHermiteGaussOrders(m(i));
                if (nx > nxMax)
                    nxMax = nx;
                end
                if (ny > nyMax)
                    nyMax = ny;
                end
            end

        end

        function labelCell = getHGOrderLabel(orderList)
            n = numel(orderList);
            labelCell = cell(1,n);

            for i = 1 : n
                labelCell{i} = sprintf('%02d', orderList(i));
            end
        end

        function str = printReadableTime(tm)
            tm = round(tm);
            SECONDS_IN_DAYS = 24*60*60;

            days = round( tm/SECONDS_IN_DAYS);
            tm = mod(tm, SECONDS_IN_DAYS);

            hours = round( tm/ (60*60));
            tm = mod(tm, 60*60);

            minutes = round(tm/60);
            seconds = mod(tm, 60);

            str = '';
            if (~days && ~hours && ~minutes && ~seconds)
                return;
            end

            if days
                str = [str num2str(days, '%d'), ' days '];
            end

            if hours
                str = [str num2str(hours, '%d'), ' hours '];
            end

            if minutes
                str = [str num2str(minutes, '%d'), ' minutes '];
            end

            if seconds
                str = [str num2str(seconds, '%d'), ' seconds '];
            end

            str = [str 'elapsed'];
        end

    end
end