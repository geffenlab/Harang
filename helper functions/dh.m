classdef dh
    %DataHelper useful functions
    
    properties
    end
    
    methods (Static)
        
        %% rasterizing
        
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
        
        %% window functions
        
        function w = make_windows(onsets, len)
            % calcualte windows around onset of trials
            % w: indices for windows
            %   num windows X len
            % onsets: onset of trial
            %   length(onsets) X 1
            % len: length of window
            w = repelem(onsets, 1, len);
            w = w + repmat(0:len-1,[length(onsets) 1]);
        end
                
        function f0_val = f0_mean(f, samp_wind)
            % calculate mean in sampling windows
            % f0_val: f0 values in the sampling windows for neurons
            % 	size(f,1) X size(window,2)
            % samp_wind: sampling window; number X window values
            f0_val = zeros(size(f,1), size(samp_wind,1));
            for i = 1 : size(samp_wind,1)
                w = samp_wind(i,:);
                f0_val(:,i) = mean(f(:,w),2);
            end
        end
        
        function f0_val = f0_min_movmean(f, samp_wind, span)
            % calculate min of sliding mean in sampling windows
            % f0_val: f0 values in the sampling windows for neurons
            % 	size(f,1) X size(window,2)
            % samp_wind: sampling window; number X window values
            % span optional (default: 22 pts; .75s for 30fps)
            % from Jia et al. Nat Protocols. 2011; 6:28-35.
            if nargin < 3
                span = 22;
            end
            f0_val = zeros(size(f,1), size(samp_wind,1));
            for i = 1 : size(samp_wind,1)
                w = samp_wind(i,:);
                f0_val(:,i) = min(movmean(f(:,w),span,2),[],2);
            end
        end
        
        function f0 = f0_apply(f, f0_val, wind)
            % apply values of f0_val to windows in f
            % f0: f0 values of size(f)
            %   f0 values outside of window (w) = f_ij
            % f: calcium trace, num neurons X time
            % f0_val: f0 values from sampled windows to 
            %   apply to the window in f, num neurons X num sampled windows
            % wind: windows to apply f0 values
            % can have 1 sampled window and apply to many windows on f0
            f0 = f;
            n_wind = size(wind,1);
            n_f0_val = size(f0_val,2);
            for i = 1 : n_wind
                if n_f0_val == 1
                    val = f0_val;
                else
                    val = f0_val(:,i);
                end
                w = wind(i,:);
                f0(:,w) = repmat(val, [1 length(w)]);
            end
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
        
        function [y, idx] = movmax(x, span, step)
            % movmax finds the max along a moving window with span
            % span default == 10
            % step default == 10
            % y: max value
            %   size(x,1) X size(x,2)/step
            % idx: index of max value in x
            %   size(x,1) X size(x,2)/step
            switch nargin
                case 1
                    span = 10;
                    step = 10;
                case 2
                    step = 10;
            end
            y = zeros(size(x,1),size(x,2)/step);
            idx = zeros(size(x,1),size(x,2)/step);
            for i = 1 : step : length(x)
                range = i : i + span - 1;
                if range(end) > length(x)
                    range = i : length(x);
                end
                ys = 1 + (i-1)/step;
                [m, mi] = max(x(:,range), [], 2);
                y(:,ys) = repmat(m, [1 length(ys)]);
                idx(:,ys) = mi + i - 1;
            end
        end
        
        function vq = change_frame_rate(v, fr0, fr)
            dur0 = size(v,2) / fr0;
            x = 1/fr0:1/fr0:dur0;
            xq = 1/fr0:1/fr:dur0;
            vq = interp1(x,v,xq);
        end
        
    end
    
end

