classdef ImRatioSeg < handle
    % radiometric intensity based segmentation:
    %   based on deviations from linear relationship between 2 fluorofore channels
    %   (from 2PM images) where I label tdtom and GCaMP6s or CsChrimson
    % Issues:
    %   1) significant tdtom signal is pickup in the green channel
    %   2) Assumes bleedthrough from the red to the green to be linear 
    %       and constant across different planes

    properties
            fname;   % file names
            fsuffix  % file suffix
            f2use;   % file to load
            sdIdx
            medMask
            verbose
            plotr
            bgate
            siz
            sig
            med_sig
            thrs_f
            thrs_std
            minpixsiz
            sd_sig
            med_siz
            oDir % output folder for fig
    end
   
   methods (Access='public')
       
        % 1) construct objects
        function obj = ImRatioSeg(fsuffix, oDir)
            
            if ~exist('fsuffix', 'var') || isempty(fsuffix) 
                obj.fsuffix = ''; 
            else
                obj.fsuffix = fsuffix; 
            end
            
            if ~exist('oDir', 'var') || isempty(oDir) 
                obj.oDir = '.'; 
            else
                obj.oDir = oDir;
                mkdir(obj.oDir);
            end
            
            if exist('images_ch2', 'dir') && exist('images', 'dir')
                obj.fname = rdir(['.', filesep, 'images_ch2', ...
                    filesep, '*', obj.fsuffix, '_02.nrrd']); 
                obj.fname = {obj.fname.name}';
            else
                obj.fname = rdir(['*', obj.fsuffix, '_Zstack.nrrd']); 
                obj.fname = {obj.fname.name}';
            end
            
            % default params
            obj.thrs_f = 25;
            obj.thrs_std = 1;
            obj.minpixsiz = 10^2;
            obj.bgate = 0;
            obj.siz = 3;
            obj.sig = 2;
            obj.sd_sig = 3;
            obj.med_siz = 3;
            obj.verbose = false;
            obj.plotr = true;
            display(obj.fname)
            
        end
        
        % 2) segment files
        function im2seg(obj, selIdx)
            
            if ~exist('selIdx', 'var') || isempty(selIdx)
               selIdx = 1:numel(obj.fname);
            end
            
            for i = selIdx
                
                fprintf(['Running file #', num2str(i), '\n'])
                
                try
                    iparams = struct(obj);
                end
                
                if contains(obj.fname{i}, '_Zstack')
                    Data = nrrdread(obj.fname{i});
                    x_ = Data(:, :, 2:2:end);
                    y_ = Data(:, :, 1:2:end);
                    iparams.fname = strrep(obj.fname{i}, '_Zstack.nrrd', '');
                else
                    x_ = nrrdread(obj.fname{i});
                    y_ = nrrdread(strrep(strrep(obj.fname{i}, '_ch2', ''), ...
                        '_02.nrrd', '_01.nrrd'));
                    iparams.fname = strrep(obj.fname{i}, '.nrrd', '');
                end
                
                [obj.medMask{i}, obj.sdIdx{i}] = ...
                    ImRatioSegment(x_, y_, iparams);
                clear Data iparams;
                
            end
            
        end
                
        % load floating images, apply median mask and save as nrrd files
        function saveresults(obj)
            
            mkdir('images_ch3')
            
            for i = 1:numel(obj.fname)
                
                fprintf(['Running file #', num2str(i), '\n'])
                
                if contains(obj.fname{i}, '_Zstack')
                    [Data, meta] = nrrdread(obj.fname{i});
                    Data = Data(:, :, 2:2:end);
                else
                    [Data, meta] = nrrdread(obj.fname{i});
                end
                
                % apply mask
                Data = double(Data).*double(obj.medMask{i});
                if contains(obj.fname{i}, '_Zstack')
                    [filename, floatFol] = split_path(obj.fname{i});
                    filename = [floatFol, strrep(filename, ...
                        '_Zstack.nrrd', '_Zstack_m.nrrd')];
                else
                    filename = strrep(strrep(obj.fname{i}, ...
                        '_02.nrrd', '_03.nrrd'), 'images_ch2', 'images_ch3');
                end
                
                % save nrrd
                nrrdWriter(filename, mat2uint16(Data, 0), ...
                    nrrdread_res(meta), [0 0 0], 'raw');
                clear Data filename;
                
            end
            
        end
        
        % save ImRatioSeg object
        function saveobj(obj)
            
            cDir = pwd;
            [filename, ~] = split_path(cDir);
            save(['imratioseg_', filename, '.mat'], 'obj');
            
        end
        
        % redo median mask for selected files (with different settings)
        function medmask_redo(obj, selIdx, med_siz, ...
                thrs_std, minpixsiz_1, minpixsiz_2)
            
            for i = selIdx
                % apply new threshold
                expr = obj.sdIdx{i};
                
                % delete small rois
                [expr_m] = split_im_into_cc(expr > 0, minpixsiz_1);
                
                expr = expr.*expr_m;
                clear expr_m;
                expr = expr > thrs_std;
                
                % median filtering
                obj.medMask{i} = medfilt_int(obj, expr, ...
                    med_siz, minpixsiz_2);
                
            end
            
        end
        
        % perform median filtering of input images
        function med_mask = medfilt_int(obj, Im, med_siz, minpixsiz)
            
            med_mask = zeros(size(Im));
            
            if numel(size(Im)) == 3
                med_mask = medfilt3(Im, [med_siz med_siz med_siz]);
            else
                med_mask = medfilt2(Im, [med_siz med_siz]);
            end
            
            % delete small rois
            [med_mask] = split_im_into_cc(med_mask > 0, minpixsiz);
            
        end
        
        % plot video of median SD mask of selected files
        function plot_slice(obj, selIdx)
            
            pi.range = [0 6];
            pi.sizY = size(obj.medMask{selIdx});
            pi.Y1 = obj.medMask{selIdx}; 
            slice3Dmatrix(obj.sdIdx{selIdx}, pi)
            
        end
        
        % plot density of defined x,y points
        function plot_density(obj, x_, y_)
            
            % density plot gate
            wpar.dgate = 0;
            % density resolution
            wpar.width = 1000;
            wpar.height = 1000;
            % smoothing params
            wpar.sig = [2 2];
            wpar.siz = [6 6];
            
            densityplot(double(x_(:))', double(y_(:))', ...
                wpar.width, wpar.height, [], wpar.sig, wpar.siz)
            
        end
        
        % plot maximun intensity projection of pairs x,y 
        function plot_MIP(obj, x_, y_)
            
            figure()
            subplot(1, 2, 1);
            imagesc(max(x_, [], 3));
            colorbar;
            title('Green Channel')
            subplot(1, 2, 2);
            imagesc(max(y_, [], 3));
            colorbar;
            title('Red Channel')
            
        end
        
        % plot video of mask and two channels
        function [Y, centers_] = plot_overlay_vid(obj, ...
            selidx, ch1_range, ch2_range, ztoignore)

            [cha_2, ~] = nrrdread(obj.fname{selidx});
            [cha_1, ~] = nrrdread(strrep(strrep(obj.fname{selidx}, ...
                'images_ch2', 'images'), '_02.', '_01.'));

            med_mask = obj.medMask{selidx};
            med_mask(:, :, 1:ztoignore) = 0;
            
            fprintf(['Extracting and expanding', ...
                ' labels (connected components)\n'])
            tic
            [med_mask] = split_im_into_cc(med_mask, ...
                obj.minpixsiz, 0);
            [Y, centers_] = labeledim_to_imperlabel(med_mask);
            toc

            % max_int_proj(med_mask)
            fprintf(['Number of labels ', num2str(size(centers_, 1)), '\n'])

            % save video
            pi.cmap = buildcmap('wr');
            pi.Y2cmap = buildcmap('wk');
            pi.Y2range = ch1_range;
            pi.Y2 = double(cha_1);
            pi.range = ch2_range;
            pi.sizY = size(cha_2);
            pi.Y1 = Y;
            pi.Y1cmap = winter(size(centers_, 1));
            pi.Y3 = centers_;
            pi.Y3(:, 3) = round(pi.Y3(:, 3));
            pi.Y3text = cf(@(x) num2str(x), ...
                    chunk2cell(1:size(centers_, 1), 1));
            pi.Y3fontsiz = 20;
            pi.vgate = 1;
            pi.vname = strrep(split_path(obj.fname{selidx}), '.nrrd', '_label');

            slice3Dmatrix(double(cha_2), pi)

        end

        % plot video of mask and two channels
        function edit_cc_and_plot_overlay_vid(obj, ...
            selidx, ch1_range, ch2_range, Y, ztoignore, labels2merge)

            [cha_2, meta] = nrrdread(obj.fname{selidx});
            [cha_1, ~] = nrrdread(strrep(strrep(obj.fname{selidx}, ...
                'images_ch2', 'images'), '_02.', '_01.'));

            if isempty(Y)
                med_mask = obj.medMask{selidx};
                med_mask(:, :, 1:ztoignore) = 0;

                fprintf(['Extracting and expanding', ...
                    ' labels (connected components)\n'])
                tic
                [med_mask] = split_im_into_cc(med_mask, ...
                    obj.minpixsiz, 0);
                [Y, centers_] = labeledim_to_imperlabel(med_mask);
                toc
            end

            % edit labels
            pre_Y = sparse(size(Y, 1), size(labels2merge, 2));
            for i = 1:numel(labels2merge)
                pre_Y(:, i) = max(Y(:, labels2merge{i}), [], 2);
            end
            Y = pre_Y;
            clear pre_Y
            
            %max_int_proj(Y)

            if ~exist('centers_', 'var') || isempty(centers_)

                % get original size
                d1 = size(cha_2, 1);
                d2 = size(cha_2, 2);
                d3 = size(cha_2, 3);

                centers_ = com(Y, d1, d2, d3);

            end

            % Plot images with these new masks
            fprintf(['Number of labels ', num2str(size(centers_, 1)), '\n'])

            % save video
            pi.cmap = buildcmap('wr');
            pi.Y2cmap = buildcmap('wk');
            pi.Y2range = ch1_range;
            pi.Y2 = double(cha_1);
            pi.range = ch2_range;
            pi.sizY = size(cha_2);
            pi.Y1 = Y;
            pi.Y1cmap = winter(size(centers_, 1));
            pi.Y3 = centers_;
            pi.Y3(:, 3) = round(pi.Y3(:, 3));
            pi.Y3text = cf(@(x) num2str(x), ...
                    chunk2cell(1:size(centers_, 1), 1));
            pi.Y3fontsiz = 20;
            pi.vgate = 1;
            pi.vname = strrep(split_path(obj.fname{selidx}), '.nrrd', '_label_edit');

            slice3Dmatrix(double(cha_2), pi)

            % collapse second dimension of labels to just one with
            %   different values per label
            Y = imperlabel_to_labeledim(Y, size(cha_2));

            % Save nrrd with these new labels
            Y = double(Y);

            if contains(obj.fname{selidx}, '_Zstack')
                [filename, floatFol] = split_path(obj.fname{selidx});
                filename = [floatFol, strrep(filename, ...
                    '_Zstack.nrrd', '_Zstack_m.nrrd')];
            else
                filename = strrep(strrep(obj.fname{selidx}, ...
                    '_02.nrrd', '_04.nrrd'), 'images_ch2', 'images_ch4');
                [~, floatFol] = split_path(filename);
                mkdir(floatFol)
            end

            % save nrrd for current and mirror image
            nrrdWriter(filename, mat2uint8(Y, 0), ...
                nrrdread_res(meta), [0 0 0], 'gzip');

            if contains(filename, '_w_')
                nrrdWriter(strrep(filename, '_w_', '_wm_'), ...
                    mat2uint8(flip(Y, 2), 0), ...
                    nrrdread_res(meta), [0 0 0], 'gzip');
            elseif contains(filename, '_wm_')
                nrrdWriter(strrep(filename, '_wm_', '_w_'), ...
                    mat2uint8(flip(Y, 2), 0), ...
                    nrrdread_res(meta), [0 0 0], 'gzip');
            end

        end

   end
   
end
