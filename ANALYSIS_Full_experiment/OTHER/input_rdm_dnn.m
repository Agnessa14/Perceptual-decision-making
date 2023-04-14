function input_rdm_dnn(corr_type,model_type)
%INPUT_RDM_DNN Compute Input RDM for the DNN from activations.
%Input: corr_type ('pearson' or 'euclidean'), model_type ('bl','b_d' or 'b') 
%Output: N num_samples x num_samples RDMs, where N=number of total DNN layers
%used (num layers x num timesteps)
%

%% load activations, define num layers and timesteps
if strcmp(model_type,'bl')
    activ_dir='/scratch/agnek95/PDM/DATA/RNN_ACTIVATIONS/activations';
    save_dir_model='02.11_2_rnn';
    num_layers=7;
    num_timesteps=8;
elseif strcmp(model_type,'b_d')
    activ_dir='/scratch/agnek95/PDM/DATA/B_D_ACTIVATIONS/activations';
    save_dir_model='B_D_net';
    num_layers=14;
    num_timesteps=1;
elseif strcmp(model_type,'b')
    activ_dir='/scratch/agnek95/PDM/DATA/B_ACTIVATIONS/activations';
    save_dir_model='B_net';
    num_layers=7;
    num_timesteps=1;
end

%% compute RDM
num_samples=60;
input_rdm=NaN(num_samples,num_samples);

for layer=1:num_layers
    l=layer-1;
    for timestep=1:num_timesteps
        t=timestep-1;
        if strcmp(model_type,'bl')
            filename_l=sprintf('ReLU_Layer_%d_Time_%d',l,t);
        elseif strcmp(model_type,'b_d') || strcmp(model_type,'b')
            filename_l=sprintf('ReLU_Layer_%d',l);
        end
        activ_path=fullfile(activ_dir,sprintf('%s_activations.mat',filename_l));    
        load(activ_path,'data');
        fprintf('Loaded activations RDM layer %d timestep %d\n',layer,timestep);
        for sample1=1:num_samples
            for sample2=1:num_samples
                if strcmp(corr_type,'pearson')
                    input_rdm(sample1,sample2)=1-corr(squeeze(data(sample1,:))',squeeze(data(sample2,:))','type','Pearson');
                elseif strcmp(corr_type,'euclidean')
                    input_rdm(sample1,sample2)= pdist2(squeeze(data(sample1,:)),squeeze(data(sample2,:)));
                end
            end
        end
    fprintf('Finished RDM layer %d timestep %d\n',layer,timestep);   
   
    save_dir=fullfile('/home/agnek95/SMST/PDM_FULL_EXPERIMENT/RESULTS_AVG/',save_dir_model,sprintf('Input_RDM_%s',corr_type));
    save(fullfile(save_dir,sprintf('%s_Input_RDM_%s.mat',filename_l,corr_type)),'input_rdm');
    fprintf('Saved RDM layer %d timestep %d\n',layer,timestep);
    end
end
end