function probe = buzsaki64probeLoc

    xShank1 = [-20.5 -16.5 -12.5 -8.5 0 8.5 12.5 16.5]'; 
    yShank1 = [140 100 60 20 0 40 80 120]';
    
    probe.chanNo  = [1:64]';
    probe.shankNo = [  ones(8,1); ...
                     2*ones(8,1); ...
                     3*ones(8,1); ...
                     4*ones(8,1); ...
                     5*ones(8,1); ...
                     6*ones(8,1); ...
                     7*ones(8,1); ...
                     8*ones(8,1)]; 
                 
    probe.x = [xShank1;
               xShank1 + 200;
               xShank1 + 400;
               xShank1 + 600;
               xShank1 + 800;
               xShank1 + 1000;
               xShank1 + 1200;
               xShank1 + 1400];
           
    probe.y = [yShank1;
               yShank1;
               yShank1;
               yShank1;
               yShank1;
               yShank1;
               yShank1;
               yShank1];

end