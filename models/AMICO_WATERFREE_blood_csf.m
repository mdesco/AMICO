classdef AMICO_WATERFREE

properties
    id, name                % id and name of the model
    dPar                    % parallel diffusivity of the tensors [units of mm^2/s]
    dPer                    % perpendicular diffusivities of the tensors [units of mm^2/s]
    dIso                    % isotropic diffusivities [units of mm^2/s]
    OUTPUT_names            % suffix of the output maps
    OUTPUT_descriptions     % description of the output maps
end


methods

    % =================================
    % Setup the parameters of the model
    % =================================
	function obj = AMICO_WATERFREE()
        global CONFIG

        % set the parameters of the model
        obj.id        = 'WATERFREE';
        obj.name      = 'Water free';
        obj.dPar      = 1.7 * 1E-3;
        obj.dIso      = [2.0 3.0] * 1E-3;
        obj.dPer      = linspace(0.1,1.0,10) * 1E-3;
        
        obj.OUTPUT_names        = { 'ICVF', 'FW' };
        obj.OUTPUT_descriptions = {'Intra-cellular volume fraction', ...
                'Isotropic free-water volume fraction'};

        % set the parameters to fit it
        CONFIG.OPTIMIZATION.SPAMS_param.mode    = 2;
        CONFIG.OPTIMIZATION.SPAMS_param.pos     = true;
        CONFIG.OPTIMIZATION.SPAMS_param.lambda  = 0;    % l1 regularization
        CONFIG.OPTIMIZATION.SPAMS_param.lambda2 = 1e-3; % l2 regularization
    end


    % ==================================================================
    % Generate high-resolution kernels and rotate them in harmonic space
    % ==================================================================
    function GenerateKernels( obj, ATOMS_path, schemeHR, AUX, idx_IN, idx_OUT )
        global CONFIG AMICO_data_path

        % Tensor compartments
        % ===================
        idx = 1;
        for i = 1:numel(obj.dPer)
            TIME = tic();
            fprintf( '\t\t- A_%03d... ', idx );

            % generate
            D = diag( [obj.dPer(i) obj.dPer(i) obj.dPar] );
            signal = obj.TensorSignal( D, schemeHR.camino );

            % rotate and save
            lm = AMICO_RotateKernel( signal, AUX, idx_IN, idx_OUT, false );
            save( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), '-v6', 'lm' )
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end


        % Isotropic compartments
        % ======================
        for i = 1:numel(obj.dIso)
            TIME = tic();
            fprintf( '\t\t- A_%03d... ', idx );

            % generate
            D = diag( [obj.dIso(i) obj.dIso(i) obj.dIso(i)] );
            signal = obj.TensorSignal( D, schemeHR.camino );

            % resample and save
            lm = AMICO_RotateKernel( signal, AUX, idx_IN, idx_OUT, true );
            save( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), '-v6', 'lm' )
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end
 
    end


    % ==============================================
    % Project kernels from harmonic to subject space
    % ==============================================
    function ResampleKernels( obj, ATOMS_path, idx_OUT, Ylm_OUT )
        global CONFIG AMICO_data_path KERNELS

        % Setup the KERNELS structure
        % ===========================
        n1 = numel(obj.dPer);
        n2 = numel(obj.dIso);
        KERNELS = {};
        KERNELS.nS       = CONFIG.scheme.nS;
        KERNELS.nA       = n1 + n2; % number of atoms
        KERNELS.A1       = zeros( [KERNELS.nS n1 181 181], 'single' );
        KERNELS.A2       = zeros( [KERNELS.nS n2], 'single' );


        % Tensors
        % =======
        idx = 1;
        for i = 1:n1
            TIME = tic();
            fprintf( '\t- A_%03d...  ', idx );

            load( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), 'lm' );
            KERNELS.A1(:,i,:,:) = AMICO_ResampleKernel( lm, idx_OUT, Ylm_OUT, false );
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end


        % Isotropic
        % =========
        for i = 1:n2
            TIME = tic();
            fprintf( '\t- A_%03d...  ', idx );

            load( fullfile( ATOMS_path, sprintf('A_%03d.mat',idx) ), 'lm' );
            KERNELS.A2(:,i,:,:) = AMICO_ResampleKernel( lm, idx_OUT, Ylm_OUT, true );
            idx = idx + 1;

            fprintf( '[%.1f seconds]\n', toc(TIME) );
        end
 
    end


    % ===========================
    % Fit the model to each voxel
    % ===========================
    function [DIRs, MAPs] = Fit( obj )
        global CONFIG
        global niiSIGNAL niiMASK
        global KERNELS bMATRIX
    
        % setup the output files
        MAPs         = zeros( [CONFIG.dim(1:3) numel(obj.OUTPUT_names)], 'single' );
        DIRs         = zeros( [CONFIG.dim(1:3) 3], 'single' );
        niiY = niiSIGNAL;
        niiY2 = niiSIGNAL;
        niiY3 = niiSIGNAL;
        niiY4 = niiSIGNAL;
        niiY5 = niiSIGNAL;
        niiY6 = niiSIGNAL;
        niiY7 = niiSIGNAL;

        residual_map = zeros( CONFIG.dim(1:3) );
        
        % number of unknown single-fiber atoms to estimate
        n1 = numel(obj.dPer);
        % number of isotropic compartments to estimate
        n2 = numel(obj.dIso);
        
        fprintf( '\n-> Fitting "%s" model to data:\n', obj.name );
        TIME = tic();
        for iz = 1:niiSIGNAL.hdr.dime.dim(4)
        for iy = 1:niiSIGNAL.hdr.dime.dim(3)
        for ix = 1:niiSIGNAL.hdr.dime.dim(2)
            if niiMASK.img(ix,iy,iz)==0, continue, end

            % read the signal
            b0 = mean( squeeze( niiSIGNAL.img(ix,iy,iz,CONFIG.scheme.b0_idx) ) );
            if ( b0 < 1e-3 ), continue, end
            y = double( squeeze( niiSIGNAL.img(ix,iy,iz,:) ) ./ ( b0 + eps ) );
            y( y < 0 ) = 0; % [NOTE] this should not happen!

            % y is the normalized signal, ndir x 1

            % find the MAIN DIFFUSION DIRECTION using DTI
            [ ~, ~, V ] = AMICO_FitTensor( y, bMATRIX );
            Vt = V(:,1);
            if ( Vt(2)<0 ), Vt = -Vt; end

            % V is 3x1 vector representing main e-vector
            
            % build the DICTIONARY
            [ i1, i2 ] = AMICO_Dir2idx( Vt );
            A = double( [ KERNELS.A1(CONFIG.scheme.dwi_idx,:,i1,i2) KERNELS.A2(CONFIG.scheme.dwi_idx,:) ] );

            % [A] contains in teh first blocs, the single fiber atoms
            % in the second block, the isotropic atoms
            
            % fit AMICO
            y = y(CONFIG.scheme.dwi_idx);
            yy = [ 1 ; y ];
            AA = [ ones(1,size(A,2)) ; A ];

            % yy contains only the normalized signal and a dummy b0=1 image
            % hence, adding a row of 1's at the beginning of the dico
            
            % estimate CSF partial volume and remove it
            % x are all the unknowns
            x = full( mexLasso( yy, AA, CONFIG.OPTIMIZATION.SPAMS_param ) );
                        
            % STORE results	
            DIRs(ix,iy,iz,:) = Vt; % fiber direction

            % single fiber contributions. the first n1 components of x
            MAPs(ix,iy,iz,1) = sum( x(1:n1) ) / ( sum(x) + eps ); % intracellular volume fraction

            % this is supposed to be the volume fraction
            % sum( x(n1+1:end) ) / (sum(x) + eps) == 1.0 - MAPs(ix, iy, iz, 1)
            
            MAPs(ix,iy,iz,2) = 1.0 - MAPs(ix,iy,iz,1); % isotropic volume fraction
            MAPs(ix,iy,iz,3) = sum( x(n1+1) ) / ( sum(x) + eps );
            MAPs(ix,iy,iz,4) = sum( x(n1+2:end) ) / ( sum(x) + eps );

            estimated_signal = AA(2:end,:)*x*b0;
            original_signal = y*b0;
                        
            xx = x;
            sum(xx);
            
            % fiber part off
            xx(1:n1) = 0;
            free_water = AA(2:end,:)*xx*b0;
            tmp = xx(n1+1);
            xx(n1+1) = 0;
            free_water_csf = AA(2:end,:)*xx*b0;
            xx(n1+1) = tmp;
            xx(n1+2) = 0;
            free_water_blood = AA(2:end,:)*xx*b0;
            
            %if x(n1+1:end) ~= 0 
            residual = norm(original_signal(:)-estimated_signal(:))/norm(original_signal(:));
            MAPs(ix,iy,iz,5) = residual;
            
            verbose = 0;
            if verbose                              
                if sum(x(n1:1:end)) >= 0.2
                    sum(x(n1:1:end))
                    original_signal'
                    free_water'
                    
                    fprintf('voxel (%d, %d, %d) - %f %f\n', ix, iy, iz, x(n1+1:end))
                    fprintf('residual - %f\n', residual);            
                    fprintf('fw - %f\n', free_water(1));
                    disp(x')
                end
            end
            
            %end
            residual = original_signal(:)-estimated_signal(:);
            % erase the coefficients in charge of isotropic compartments
            tmp = x(n1+1);
            % erase the FW blood component
            x(n1+1) = 0;
            FWblood_estimated_signal = AA(2:end,:)*x*b0;
            % put back the blood component and erase the CSF
            x(n1+1) = tmp;
            x(n1+2) = 0;
            FWcsf_estimated_signal = AA(2:end,:)*x*b0;
            % erase both
            x(n1+1:end) = 0;            
            FiberFraction_estimated_only_signal = AA(2:end,:)*x*b0 + residual;
            FiberFraction_estimated_signal = AA(2:end,:)*x*b0;
            % re-generate the free-water corrected DWI signal
            %signal_fw_corrected = original_signal - free_water;
            %signal_blood_corrected = original_signal - free_water_blood;
            %signal_csf_corrected = original_signal - free_water_csf;

            niiY.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = estimated_signal;
            niiY2.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = FiberFraction_estimated_only_signal;
            niiY3.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = FiberFraction_estimated_signal;;
            niiY4.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = FWblood_estimated_signal;
            niiY5.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = FWcsf_estimated_signal;
            niiY6.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = FWblood_estimated_signal+residual;
            niiY7.img(ix,iy,iz,CONFIG.scheme.dwi_idx) = FWcsf_estimated_signal+residual;
        end
        end
        end
        save_untouch_nii(niiY, fullfile(CONFIG.OUTPUT_path,'dwi_estimated.nii'));
        save_untouch_nii(niiY2, fullfile(CONFIG.OUTPUT_path,'dwi_FiberComparment_only.nii'));
        save_untouch_nii(niiY3, fullfile(CONFIG.OUTPUT_path,'dwi_FiberComparment_est.nii'));
        save_untouch_nii(niiY4, fullfile(CONFIG.OUTPUT_path,'dwi_fw_blood_corrected.nii'));
        save_untouch_nii(niiY5, fullfile(CONFIG.OUTPUT_path,'dwi_fw_csf_corrected.nii'));
        save_untouch_nii(niiY6, fullfile(CONFIG.OUTPUT_path,'dwi_fw_blood_corrected_r.nii'));
        save_untouch_nii(niiY7, fullfile(CONFIG.OUTPUT_path,'dwi_fw_csf_corrected_r.nii'));

        
        
        TIME = toc(TIME);
        fprintf( '   [ %.0fh %.0fm %.0fs ]\n', floor(TIME/3600), floor(mod(TIME/60,60)), mod(TIME,60) )

    end


    % ================================================================
    % Simulate signal according to tensor model (1 fiber along z-axis)
    % ================================================================
    function [ signal ] = TensorSignal( obj, D, XYZB )
        nDIR   = size( XYZB, 1 );
        signal = zeros( nDIR, 1 );
        for d = 1:nDIR
            signal(d) = exp(-XYZB(d,4) * XYZB(d,1:3) * D * XYZB(d,1:3)');
        end
    end

end

end
