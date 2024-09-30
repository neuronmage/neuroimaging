clear all
for sub = [01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18]  
    for set = [1 2 3]
        for run = [1 2]
            for fract = [1 2 3 4 5 6 7 8]
                %Turn the current sub/set/run 
                curr_fract_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/RL2_0%d/set%d_run%d_fract%d_corr.txt', sub, set, run, fract);
                curr_fract_file
                curr_stim_data = csvread(curr_fract_file);
                curr_stim_data = curr_stim_data';
                curr_save_file_name = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/RL2_0%d/set%d_run%d_memstr_fract%d', sub, set, run, fract);
                runanalysis(curr_stim_data, 1, 0.25, curr_save_file_name, 0.005, 0);

                load(curr_save_file_name)

                curr_fract_pmode_save_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/RL2_0%d/set%d_run%d_pmode_fract%d.txt', sub, set, run, fract);
                save(curr_fract_pmode_save_file,'pmode','-ASCII')

                curr_fract_p05_save_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/RL2_0%d/set%d_run%d_p05_fract%d.txt', sub, set, run, fract);
                save(curr_fract_p05_save_file,'p05','-ASCII')

                curr_fract_p95_save_file = sprintf('/home/data/madlab/data/mri/rev_lrn/pilot_data/behavioral/learning_analysis/RL2_0%d/set%d_run%d_p95_fract%d.txt', sub, set, run, fract);
                save(curr_fract_p95_save_file,'p95','-ASCII')
            end    
        end
    end
end

