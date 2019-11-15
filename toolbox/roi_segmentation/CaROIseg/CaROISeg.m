classdef CaROISeg < handle
   
    % Methods info:
    % 1) CaROISeg: initialize 'obj' and loads default 'CNMFSetParms'
    % 2) updateParams: updates fields within 'CNMFSetParms'
    % 3) loaddata: load memory map file
    % 4) makepatch: split volume/planes into XYZ patches to run in paralel
    % 5) initComponents_patches: data preprocessing and initialization patches
    % 6) initComponents_part1: generate RESULTS for each segment separatelly
    % 7) initComponents_part2: compile RESULTS and finish segmentation
    % 8) initComponents_part3: compile RESULTS and finish segmentation
    % 9) refineComponents: manual refinement
    % 10) updateSpatial: update spatial components
    % 11) updateTemporal: update temporal components
    % 12) updatecenter: update centers
    % 13) merge: merge found components
    % 14) lcgen: compute local corr on patches
    % 15) thresholdA: threshold components
    % 16) residual: compute the residual whole area/volume
    % 17) bgsub: compute background substraction whole area/volume
    % 18) medsub: compute median substraction whole area/volume
    % 19) rawsignal: compute additional traces
    % 20) rawsignal_red: compute additional traces (red only)
    % 21) bas_estimate: estimate baseline activity
    % 22) snapshot:take the snapshot of current results
    % 23) populate_obj: Uses an ROI snapshot to populate obj before using the object
    % 24) orderROIs: order_ROIs
    % 25) pca_Yr: PCA related    

   properties
        data;       % map data file
        data_;      % map data file (extra channel)
        sizeY;      % Y original size
        patches;    % patches generated
        A;          % spatial components of neurons 
        C;          % temporal components of neurons 
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        S;          % spike counts 
        options;    % options for model fitting
        P;          % some estimated parameters
        center;     % center of ROIs
        RESULTS;    % Results of the CNMF algorithm on individual patches
        Yr;         % Residual whole data
        lc;         % localcorr map 2D
        % ****************************   signalExtraction_int   *****************************************
        inferred;   % spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
        filtered;   % spatial weighting on the ROI pixels, unmixing, background substraction
        raw;        % uniform spatial weighting on the ROI pixels, background substraction
        filtered_;  % spatial weighting on the ROI pixels, background substraction on extra channel
        raw_;       % uniform spatial weighting on the ROI pixels, background substraction on extra channel
        % ****************************      getMode4FxT       ********************************************
        bas;        % mode of activity
        sn;         % estimated sd (from gaussian fit)
        bas_p;      % estimated mode (from gaussian fit)      
   end
   
   methods (Access = 'public')
       
        % 1) constructor and options setting
        function obj = CaROISeg(varargin)
            
            obj.options = CNMFSetParms(); 
            obj.P = struct('p', 2); 
            if nargin>0
                obj.options = CNMFSetParms(obj.options, ...
                    varargin{:}); 
            end
            
        end
        
        % 2) update parameters
        function updateParams(obj, varargin)
            
        	obj.options = CNMFSetParms(obj.options, ...
                varargin{:});
            
        end
        
        % 3) memmap mat file
        function loaddata(obj, filename)
           
            if ~contains(filename, '.mat')
                filename = [filename, '.mat'];
            end
            
            if ~exist(filename, 'file')
                % generate the memmap mat file
                fprintf('Error file does not exist\n')
            else
                % memmap *.mat file
                obj.data = matfile(filename, 'Writable', true);
            end
            
            obj.sizeY = obj.data.sizY;
            obj.options = CNMFSetParms(obj.options, ...
                'd1', obj.sizeY(1), 'd2', obj.sizeY(2));
            
            % update third dimension size
            if length(obj.sizeY) == 4
                obj.options = CNMFSetParms(obj.options, ...
                    'd3', obj.sizeY(3));
            else
                obj.options = CNMFSetParms(obj.options, ...
                    'd3', 1);
            end
            
        end
        
        % 4) construct patches
        function makepatch(obj, patch_size, ...
                overlap, patchtype, seg_idx)
            
            if ~exist('patchtype', 'var') || ...
                    isempty(patchtype)
                patchtype = 0;
            end
            
            if ~exist('seg_idx', 'var') || ...
                    isempty(seg_idx)
                seg_idx = [];
            end
            
            if patchtype == 0
                obj.patches = construct_patches(...
                    obj.sizeY(1:end-1), ...
                    patch_size, overlap);
            else
                fprintf('Using segment idx to split patches\n')
                obj.patches = construct_patches_seg(...
                    obj.sizeY(1:end-1), seg_idx);
            end
            
        end
        
        % 5) data preprocessing and initialization patches
        function initComponents_patches(obj, ...
                K, tau, p)
            
            [obj.A, obj.b, obj.C, obj.f, obj.S, obj.P] = ...
                run_CNMF_patches_int(obj.data, K, ...
                obj.patches, tau, p, obj.options);
            
            % update centers
            updatecenter(obj)
            
        end
        
        % 6) generate RESULTS for each segment separatelly
        function initComponents_part1(obj, ...
                K, tau, p, patchidx, filename)
            
            run_CNMF_patches_int_runperseg(obj.data, ...
                K, obj.patches, tau, p, obj.options, ...
                patchidx, filename);
            
        end    
        
        % 7) compile RESULTS and finish segmentation
        function initComponents_part2(obj, ...
                K, tau, p, filename)
            
            [obj.A, obj.b, obj.C, obj.f, obj.S, obj.P] = ...
                run_CNMF_patches_int_compile(obj.data, ...
                K, obj.patches, tau, p, obj.options, filename);
            
            % update centers
            updatecenter(obj)
            
        end    
        
        % 8) compile RESULTS and finish segmentation
        function initComponents_part3(obj, ...
                K, tau, p, patchidx, ...
                filename, rfisuffix)
            
            run_CNMF_patches_selected_segment(obj.data, K, ...
                obj.patches, tau, p, obj.options, patchidx, ...
                filename, rfisuffix);
            
        end  
        
        % 9) manual refinement
        function refineComponents(obj, sx)
            
            if ~exist('sx', 'var') || ...
                    isempty(sx)
                sx = 5;
            end
            
            Y = obj.data.Y;
            [obj.A, obj.C, obj.center] = ...
                manually_refine_components(...
                Y, obj.A, obj.C, obj.center, ...
                obj.lc, sx, obj.options);
            
        end
        
        % 10) update spatial components
        function updateSpatial(obj)
            
            if strcmpi(obj.options.spatial_method, 'regularized')
                A_ = [obj.A, obj.b];
            else
                A_ = obj.A;
            end
            
            [obj.A, obj.b, obj.C] = ...
                update_spatial_components(obj.data, ...
                obj.C, obj.f, A_, obj.P, obj.options);
            
            % update centers
            updatecenter(obj)
            
        end
        
        % 11) update temporal components
        function updateTemporal(obj)
            
            [obj.C, obj.f, obj.P, obj.S] = ...
                update_temporal_components(...
                obj.data, obj.A, obj.b, ...
                obj.C, obj.f, obj.P, obj.options);
            
            % update centers
            updatecenter(obj)
            
        end
        
        % 12) update centers
        function updatecenter(obj)
            
            if length(obj.sizeY) == 3
                obj.center = com(obj.A, ...
                    obj.sizeY(1), obj.sizeY(2)); 
            else
                obj.center = com(obj.A, ...
                    obj.sizeY(1), ...
                    obj.sizeY(2), obj.sizeY(3));
            end
            
        end
                       
        % 13) merge found components
        function [nr, merged_ROIs] = merge(obj)
            
            [obj.A, obj.C, nr, merged_ROIs, obj.P, obj.S] = ...
                merge_components(obj.data.Y, obj.A, ...
                obj.b, obj.C, obj.f, obj.P, obj.S, obj.options);
                        
            % update centers
            updatecenter(obj)
            
        end
        
        % 14) compute local corr on patches
        function lcgen(obj, iparams)
            
            if ~exist('iparams', 'var')
                iparams = [];
            end
            
            if length(obj.sizeY) == 3 
                obj.lc = neighcorr_2D(obj.data, iparams);
            else
                obj.lc = neighcorr_3D(obj.data, iparams);
            end
            
        end
        
        % 15) threshold components
        function thresholdA(obj)
            
            obj.A = threshold_components(obj.A, obj.options);
            
        end
        
        % 16) compute the residual whole area/volume
        function Yr = residual(obj, Yr)
            
            if length(size(Yr)) > 2
                Yr = Yr - reshape(obj.A*obj.C ...
                    + obj.b*obj.f, obj.sizeY);
            else
                Yr = Yr - obj.A*obj.C - obj.b*obj.f;
            end
            
        end
        
        % 17) compute background substraction whole area/volume
        function Yb = bgsub(obj, Yr)
            
            if length(size(Yr)) > 2
                Yb = Yr - reshape(obj.b*obj.f, obj.sizeY);
            else
                Yb = Yr - obj.b*obj.f;
            end
            
        end
        
        % 18) compute median substraction whole area/volume
        function Ymed = medsub(obj, Yr)
            Ymed = bsxfun(@minus, Yr, median(Yr, 4));
        end
        
        % 19) compute additional traces
        function rawsignal(obj, filename)
            
            if ~contains(filename, '.mat')
                filename = [filename, '.mat'];
            end
            
            % map red channel if it exist
            if exist('filename', 'var') ...
                    && ~isempty(filename)
                
                try
                    obj.data_ = matfile(filename, 'Writable', true);
                catch
                    fprintf('Red channel does not exist\n')
                end
                
            end
            
            [obj.inferred, obj.filtered, obj.raw, ...
                obj.filtered_, obj.raw_] = ...
                signalExtraction_int(obj.data, ...
                obj.data_, obj.A, obj.C, ...
                obj.b, obj.f, obj.options.d1, ...
                obj.options.d2, obj.options.d3, ...
                obj.options);
            
        end
        
        % 20) compute additional traces (red only)
        function rawsignal_red(obj, filename)
            
            if ~contains(filename, '.mat')
                filename = [filename, '.mat'];
            end
            
            % map red channel if it exist
            if exist('filename', 'var') ...
                    && ~isempty(filename)
                
                try
                    obj.data_ = matfile(filename, 'Writable', true);
                catch
                    fprintf('Red channel does not exist\n')
                end
                
            end
            
            [obj.filtered_, obj.raw_] = ...
                signalExtraction_int_2(obj.data, ...
                obj.data_, obj.A, obj.C, obj.b, ...
                obj.f, obj.options.d1, ...
                obj.options.d2, obj.options.d3, ...
                obj.options);
            
        end
        
        % 21) estimate baseline activity
        function bas_estimate(obj)
            
            % Estimate baseline activity (to be deprecated)
            [obj.bas, obj.sn, obj.bas_p] = ...
                getMode4FxT(obj.filtered.dfof);
            
        end
                
        % 22) take the snapshot of current results
        function roi = snapshot(obj)
            
            roi.A = obj.A; 
            roi.C = obj.C;
            roi.b = obj.b; 
            roi.f = obj.f;
            roi.P = obj.P;
            roi.size = obj.sizeY;
            roi.center = obj.center;
            roi.options = obj.options;
            roi.inferred = obj.inferred;
            roi.filtered = obj.filtered;
            roi.raw = obj.raw;
            roi.filtered_ = obj.filtered_;
            roi.raw_ = obj.raw_;
            roi.bas = obj.bas;
            roi.sn = obj.sn;
            roi.bas_p = obj.bas_p;
            
        end
        
        % 23) Uses an ROI snapshot to populate obj before using the object
        function populate_obj(obj, roi)
            
            obj.A = roi.A; 
            obj.C = roi.C;
            obj.b = roi.b; 
            obj.f = roi.f;
            obj.P = roi.P;
            obj.sizeY = roi.size;
            obj.center = roi.center;
            obj.options = roi.options;
            obj.inferred = roi.inferred;
            obj.filtered = roi.filtered;
            obj.raw = roi.raw;
            obj.filtered_ = roi.filtered_;
            obj.raw_ = roi.raw_;
            obj.bas = roi.bas;
            obj.sn = roi.sn;
            obj.bas_p = roi.bas_p;
            
        end
                
        % 24) order_ROIs 
        function [srt] = orderROIs(obj)
            
            [obj.A, obj.C, obj.S, obj.P, srt] = ...
                order_ROIs(obj.A, obj.C, obj.S, obj.P);
            
            % update centers
            updatecenter(obj)
            
        end
        
        % 25) PCA related
        function [rPCA] = pca_Yr(obj, pcacomp_n)
            
            % Would only work for small volumes
            %   (preferably single planes)
            Y = obj.data.Yr;
            
            if size(Y, 1) < size(Y, 2)
                Y = Y';
            end
            
            [rPCA.coeff, rPCA.score, rPCA.latent, ...
                rPCA.tsquared, rPCA.explained, rPCA.mu] = ...
                pca(Y , 'NumComponents', pcacomp_n);
            
        end
        
   end
   
end
