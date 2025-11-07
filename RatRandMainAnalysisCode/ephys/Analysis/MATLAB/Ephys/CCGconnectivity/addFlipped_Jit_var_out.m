function jit_var_out = addFlipped_Jit_var_out(jit_var_out)

    lenConstant = length(jit_var_out);

    for i = 1:lenConstant

        jit_var_out(i + lenConstant).cell_pair         = fliplr(jit_var_out(i).cell_pair);
        jit_var_out(i + lenConstant).GSPExc            = fliplr(jit_var_out(i).GSPExc);
        jit_var_out(i + lenConstant).GSPInh            = fliplr(jit_var_out(i).GSPInh);
        jit_var_out(i + lenConstant).pvalE             = fliplr(jit_var_out(i).pvalE);
        jit_var_out(i + lenConstant).pvalI             = fliplr(jit_var_out(i).pvalI);
        jit_var_out(i + lenConstant).ccgR              = fliplr(jit_var_out(i).ccgR);
        jit_var_out(i + lenConstant).tR                = fliplr(jit_var_out(i).tR);
        jit_var_out(i + lenConstant).LSPExc            = fliplr(jit_var_out(i).LSPExc);
        jit_var_out(i + lenConstant).LSPInh            = fliplr(jit_var_out(i).LSPInh);
        jit_var_out(i + lenConstant).JBSIE             = fliplr(jit_var_out(i).JBSIE);
        jit_var_out(i + lenConstant).JBSII             = fliplr(jit_var_out(i).JBSII);
        jit_var_out(i + lenConstant).cell1type         = jit_var_out(i).cell2type;
        jit_var_out(i + lenConstant).cell2type         = jit_var_out(i).cell1type;
        jit_var_out(i + lenConstant).cell1shank        = jit_var_out(i).cell2shank;
        jit_var_out(i + lenConstant).cell2shank        = jit_var_out(i).cell1shank;
        jit_var_out(i + lenConstant).epoch             = jit_var_out(i).epoch;
        jit_var_out(i + lenConstant).cell1_spike_times = jit_var_out(i).cell2_spike_times;
        jit_var_out(i + lenConstant).cell2_spike_times = jit_var_out(i).cell1_spike_times;

    end
    
end

