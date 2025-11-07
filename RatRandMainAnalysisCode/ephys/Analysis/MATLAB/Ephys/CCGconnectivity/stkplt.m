function [bias0,scale0] = stkplt(action,arg2)
%STKPLT stacks plots [bias0,scale0] = stkplt(action,arg2)
%
%
% bias   - offsets on y-axis
% scale  - re-scaling factors
%
% Author: Bert deVries, at some point in the late 90'

if nargin<1, action = 'init'; end


%%% find plots, get rid of vertical lines
hp = findobj(gca,'type','line'); hplots=[];
for i=1:length(hp)
   xd = get(hp(i),'xdata');
   if any(xd~=xd(1)), hplots=[hplots,hp(i)]; end %if
end
nsig = length(hplots); 


if strcmp(action,'play')
   h = findobj(gcf,'type','uicontrol','string','play');
   for i=1:length(h)
      if get(h(i),'value')
         set(h(i),'value',0);
         y =  get(hplots(nsig+1-i),'YData');
         soundsc(y-mean(y));
      end%if
   end%for

else

grp = 1:nsig; ngrp = max(grp);

space = 0.1/ngrp;  
range = (1-(ngrp+1)*space)/ngrp;

if isempty(getappdata(gca,'stack'))
  setappdata(gca,'stack','off');
end

if strcmp(getappdata(gca,'stack'),'off')
    
  % stack
  bias = zeros(nsig,1); scale = zeros(nsig,1);
  for i=1:nsig
    y = get(hplots(i),'YData');
    ymax = nanmax(y); ymin = nanmin(y);
    offset = grp(i)*space + (grp(i)-1)*range;
    scale(i) = range/(ymax-ymin + eps); % scale to [0->range]
    bias(i) = offset - scale(i)*ymin;
    set(hplots(i),'YData',scale(i)*y+bias(i), 'Userdata',[scale(i) bias(i)]);
  end
  setappdata(gca,'stack','on');
  set(gca,'YTick',[]);

  % play buttons
  if nargin>1, if arg2==1

    set(gcf,'Units','norm');
    axispos = get(gca,'pos');
      
    xlim = get(gca,'xlim'); 
    ylim = get(gca,'ylim'); bias(nsig+1) = ylim(2);
    for i=1:nsig
      uicontrol('Style','pushbutton','String','play',...
        'Callback','stkplt(''play'')',...
        'Units','norm',...
        'Position',[0.01,axispos(2)+(bias(i)-ylim(1))/(ylim(2)-ylim(1))*axispos(4)-.035,.09,.07]);
    end%for
  end,end%if nargin

else

  % unstack
  delete(findobj(gcf,'style','pushbutton','string','play'));
  for i=1:nsig
    stack = get(hplots(i),'Userdata');
    if ~isempty(stack)
      scalei = stack(1); biasi = stack(2);
      y =  get(hplots(i),'YData');
      set(hplots(i),'YData',(y-biasi)/scalei,'Userdata',[]);
    end
  end
  setappdata(gca,'stack','off');
  set(gca,'YTickMode','auto');
end    

if nargout, bias0 = bias; scale0 = scale; end

end


function y = nanmax(x);
%
y = max(x(~isnan(x)));

function y = nanmin(x);
%
y = min(x(~isnan(x)));


%%%
%%% left overs
%%%

% don't use plots with all NaN's
% for i=length(hplots):-1:1
%   if all(isnan(get(hplots(i),'YData')))
%     hplots(i)=[];
%   end
% end    

