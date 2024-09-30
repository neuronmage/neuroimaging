clear all
transfer = {'near', 'far'};
%sub = [01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18]
for sub = [04]  
    for run = [1 2 3 4 5 6]
        for curr_trans = 1:length(transfer)       
            curr_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/all_levels/RL2_%03d/run%d_%s_acc.txt', sub, run, transfer{curr_trans});
            curr_file
            curr_stim_data = csvread(curr_file);
            curr_stim_data = curr_stim_data';
            curr_save_file_name = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/all_levels/RL2_%03d/run%d_%s_memstr', sub, run, transfer{curr_trans});
            runanalysis(curr_stim_data, 1, 0.25, curr_save_file_name, 0.5, 0);

            load(curr_save_file_name)

            curr_pmode_save_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/all_levels/RL2_%03d/run%d_%s_pmode.txt', sub, run, transfer{curr_trans});
            save(curr_pmode_save_file,'pmode','-ASCII')

            curr_p05_save_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/all_levels/RL2_%03d/run%d_%s_p05.txt', sub, run, transfer{curr_trans});
            save(curr_p05_save_file,'p05','-ASCII')

            curr_p95_save_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/all_levels/RL2_%03d/run%d_%s_p95.txt', sub, run, transfer{curr_trans});
            save(curr_p95_save_file,'p95','-ASCII')  
        end         
    end
end
%You will need to add the runanalysis function folder to the path --
%navigate to the folder containing the function then rerun -- add to path
