classdef dh
    %DataHelper useful functions 
    
    properties
    end
    
    methods (Static)
        
        function raster = rasterize(data, dur, onsets)
            %Rasterize creates a raster from data
            %   raster is n x length(onsets) x dur
            raster = zeros(size(data,1), dur, length(onsets));
            for i = 1 : length(onsets)
                raster(:,:,i) = data(:, onsets(i):onsets(i)+dur-1);
            end
        end
        
        function trans = find_transitions(order)
            if size(order,1) == 1
                order = order';
            end
            trans = [order(1:end-1) order(2:end)];
        end
        
        function [ uT ] = get_unique_transitions(trans)
            uT = unique(trans,'rows');
            uT = sortrows(uT,[1,2]);
        end
        
        function [ raster, num_trans ] = rasterize_by_trans(data, ...
                trans, trans_u)
            %Rasterize_by_trans sorts responses by unique transitions...
            %and averages the responses across frames
            %   raster is #neurons x stimulus duration (ie size(data,2))...
            %       x #unique transitions
            %   num_trans - num transitions
            raster = zeros(size(data,1), size(data,2), length(trans_u));
            num_trans = zeros(1,length(trans_u));
            for t = 1:length(trans_u)
                rows = find(trans(:,1)==trans_u(t,1) & ...
                    trans(:,2)==trans_u(t,2));
                num_trans(t) = length(rows);
                raster(:,:,t) = mean(data(:,:,rows),3);
            end
        end
        
        function new_dat = split_concat(dat)
             % opposite of rasterize
             n = size(dat,3);
             new_dat = zeros(size(dat,1),n*size(dat,2));
             for i = 1 : n
                 on = (i-1)*size(dat,2)+1;
                 off = on+size(dat,2)-1;
                 new_dat(:,on:off) = dat(:,:,i);
             end
        end
        
        function resps = max_in_windows( responses, duration, onsets )
            %max_in_windows finds the max responses in each window
            %   window = onset : onset + duration - 1
            resps = zeros(size(responses,1),length(onsets));
            for i = 1 : size(responses,1)
                for j = 1 : length(onsets)
                    resps(i,j) = max(...
                        responses(i,onsets(j) : onsets(j)+duration-1));
                end
            end
        end
        
        %% delta F over F0 calculations
        % Jia et al. Nat Protocols. 2011; 6:28-35.
        
        function df = deltaf_f0(f, f0_wind, wind_fun, varargin)
            % calculate delta F / F0
            % f calcium trace F, neurons X time
            % f0_wind window for f0, neurons X 2,   [[w1_on w1_off];
            %                                        [w2_on w2_off];...
            % wind_fun function to run over window
            % varargin for window function
            f0 = wind_fun(f, [f0_wind(:,1) f0_wind(:,2)], varargin{:});
            df = (f - f0) ./ f0;
        end
        
        %% calculate windows
        
        function w = window_pretrial(onsets, len)
            % onsets: length(onsets) X 1; onset of trial
            % len: length of window
            % w: index of windows
            w = repelem(onsets, 1, len);
            w = w + repmat(-1*len:-1,[length(onsets) 1]);
        end
        
        %% window functions
        
        function f0 = f0_mean(f, window)
            rep_size = [1 size(f,2)];
            f0 = repmat(mean(f(:,window), 2), rep_size);
        end
        
        function f0 = f0_mode(f, window)
            rep_size = [1 size(f,2)];
            f0 = repmat(mode(f(:,window), 2), rep_size);
        end
        
        function f0 = f0_min_moving_avg(f, window, span)
            % minimum of moving-averaged values
            % span optional (default: 22; ~0.75s for 30fps)
            if nargin < 3
                span = 22;
            end
            n_neuron = size(f,1);
            f0 = zeros(n_neuron, 1);
            for i = 1 : n_neuron
                f0(i) = min(smooth(f(i, window), span));
            end
        end
        
        function f0 = f0_min_boxcar_mean(f, window, span)
            % minimum of moving-averaged values
            % span optional (default: 22)
            if nargin < 3
                span = 22;
            end
            f0 = min(movmean(f(:,window), span, 2),[],2);
        end
        
        function f0 = f0_min_smooth(f, window, span)
            % from Jia et al.
            
        end
        
        %% math
        
        function y = denoise(x, span, smooth_fun)
            % denoise using moving average with coeff = 1/span
            % x can be neurons X time
            % span optional (default: 22 pts; .75s for 30fps)
            % smooth_fun optional (default: 'moving')
            def_span = 22; def_smooth_fun = 'moving';
            switch nargin
                case 1
                    span = def_span; smooth_fun = def_smooth_fun;
                case 2
                    smooth_fun = def_smooth_fun;
            end
            y = zeros(size(x));
            n_neuron = size(x, 1);
            for i = 1 : n_neuron
                y(i,:) = smooth(x(i,:), span, smooth_fun);
            end
        end
        
    end
    
end

