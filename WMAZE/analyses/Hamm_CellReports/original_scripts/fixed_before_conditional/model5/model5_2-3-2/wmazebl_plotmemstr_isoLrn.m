clear all
trl_type = {'S1F1' 'S1F2' 'S2F1' 'S2F2' 'S1_C' 'S2_C'};
for sub = [1 2 3 4 5 6 7 8 9]
    
    for curr_trl_type = 1:length(trl_type) 
        curr_file = sprintf('/home/data/madlab/data/behav/CAT1/CAT1_%03d/%s_corr_set%d.csv',sub,trl_type{curr_trl_type},curr_set);
        curr_stim_data = csvread(curr_file);
        curr_stim_data = curr_stim_data';
        curr_save_file_name = sprintf('/home/data/madlab/data/behav/CAT1/CAT1_%03d/%s_memstr_set%d',sub,trl_type{curr_trl_type},curr_set);
        runanalysis(curr_stim_data, 1, 0.5, curr_save_file_name, 0.005, 0);

        load(curr_save_file_name)

        curr_trl_type_pmode_save_file = sprintf('/home/data/madlab/data/behav/CAT1/pilot/CAT1_%03d/MRthesis/%s_pmode_set%d.txt',sub,trl_type{curr_trl_type},curr_set);
        save(curr_trl_type_pmode_save_file,'pmode','-ASCII')

        curr_trl_type_p05_save_file = sprintf('/home/data/madlab/data/behav/CAT1/pilot/CAT1_%03d/MRthesis/%s_p05_set%d.txt',sub,trl_type{curr_trl_type},curr_set);
        save(curr_trl_type_p05_save_file,'p05','-ASCII')

        curr_trl_type_p95_save_file = sprintf('/home/data/madlab/data/behav/CAT1/pilot/CAT1_%03d/MRthesis/%s_p95_set%d.txt',sub,trl_type{curr_trl_type},curr_set);
        save(curr_trl_type_p95_save_file,'p95','-ASCII')
    end
end

