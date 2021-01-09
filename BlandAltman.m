function [means,diffs,meanDiff,CR] = BlandAltman(var1, var2, flag,i,type)
 
    %%%Plots a Bland-Altman Plot
    %%%INPUTS:
    %%% var1 and var2 - vectors of the measurements
    %%% flag - how much you want to plot
        %%% 0 = no plot
        %%% 1 = just the data
        %%% 2 = data and the difference and CR lines
    %%%
    %%%OUTPUTS:
    %%% means = the means of the data
    %%% diffs = the raw differences
    %%% meanDiff = the mean difference
    %%% CR = the 2SD confidence limits
    %%% linfit = the paramters for the linear fit
    
    means = mean([var1;var2]);
    diffs = var1-var2;
    
    meanDiff = mean(diffs);
    sdDiff = std(diffs);
    CR = [meanDiff + 1.96 * sdDiff, meanDiff - 1.96 * sdDiff]; %%95% confidence range
    
    %%% plot results unless flag is 0
    if flag ~= 0
        plot(means,diffs,'o')
        hold on
        if flag > 1
            plot(means, ones(1,length(means)).*CR(1),means, ones(1,length(means)).*CR(2),'r-'); %%%plot the upper CR
            %plot(means, ones(1,length(means)).*CR(2),'Color',[0.8500 0.3250 0.0980],'LineStyle','-'); %%%plot the lower CR
            plot(means, ones(1,length(means)).*meanDiff,'k'); %%% plot mean
        end
        title("Bland-Altman Plot ID #"+num2str(i),type)
        ylabel("Difference between filtered and original signal")
        xlabel("Mean of filtered and original signal")
        legend('Differences', strcat('Upper Confidence Range = ',num2str(CR(1)),'s'),strcat('Lower Confidence Range = ',num2str(CR(2)),'s'),strcat('Mean = ',num2str(meanDiff),'s'), 'location','eastoutside')
    end
    
end