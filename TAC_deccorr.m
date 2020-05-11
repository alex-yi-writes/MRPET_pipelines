close all;
Appointments = [1 2 1 2 1 1 2];
t0_frame_bsl    = 95;
t0_frame_task   = 115;

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_InFlow.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Baseline.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ paths.TACs num2str(IDs(id)) num2str(d) '_TACDATA_Task.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            
            Lengths=[10*ones(30,1); 60*ones(55,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=300*ones(11,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            try
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; (TACDATA_Baseline.CerC.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.CerC.Bilateral.tac).*(2^(t0_frame_task/109))];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; (TACDATA_Baseline.Put.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Put.Bilateral.tac).*(2^(t0_frame_task/109))];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; (TACDATA_Baseline.Caud.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Caud.Bilateral.tac).*(2^(t0_frame_task/109))];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                ylim([0 45000]);
                clear titlestring
                if Appointments(id) == 1
                    titlestring = ['ID ' num2str(IDs(id)) ', morning appointment, Session type ' num2str(d)];
                else
                    titlestring = ['ID ' num2str(IDs(id)) ', afternoon appointment, Session type ' num2str(d)];
                end
                title(titlestring,'Fontsize',20,'Fontweight','bold')
                print('-dpdf','-bestfit', fullfile(paths.figures, [ num2str(IDs(id)) num2str(d) '_TAC_dc.pdf']));
                
            catch
                Cer=[vertcat(nan(9,1),TACDATA_InFlow.CerC.Bilateral.tac); (TACDATA_Baseline.CerC.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.CerC.Bilateral.tac).*(2^(t0_frame_task/109))];
                Put=[vertcat(nan(9,1),TACDATA_InFlow.Put.Bilateral.tac); (TACDATA_Baseline.Put.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Put.Bilateral.tac).*(2^(t0_frame_task/109))];
                Caud=[vertcat(nan(9,1),TACDATA_InFlow.Caud.Bilateral.tac); (TACDATA_Baseline.Caud.Bilateral.tac).*(2^(t0_frame_bsl/109)); (TACDATA_Task.Caud.Bilateral.tac).*(2^(t0_frame_task/109))];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                ax = gca; ax.YAxis.Exponent = 0;
                ylim([0 45000]);
                clear titlestring
                if Appointments(id) == 1
                    titlestring = ['ID ' num2str(IDs(id)) ', morning appointment, Session type ' num2str(d)];
                else
                    titlestring = ['ID ' num2str(IDs(id)) ', afternoon appointment, Session type ' num2str(d)];
                end
                title(titlestring,'Fontsize',20,'Fontweight','bold')
                print('-dpdf','-bestfit', fullfile(paths.figures, [ num2str(IDs(id)) num2str(d) '_TAC_dc.pdf']));
                
            end

        end
    end
end
