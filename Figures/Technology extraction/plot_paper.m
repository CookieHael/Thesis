h_plot=get(gcf,'children');
for kt=1:length(h_plot),
  h_objects=get(h_plot(kt),'children');
  if strcmp(get(h_plot(kt),'Type'),'axes')
    set(h_plot(kt),'fontsize',20);
    set(h_plot(kt),'fontweight','bold');
    set(h_plot(kt),'linewidth',2);
    h_axes=get(h_plot(kt),'children');
    for it=1:length(h_axes),
      if strcmp(get(h_axes(it),'Type'),'line'),
        set(h_axes(it),'LineWidth',3);
      end
      if strcmp(get(h_axes(it),'Type'),'text')
        set(h_axes(it),'fontsize',20);
        set(h_axes(it),'fontweight','bold');
      end
    end
    h_label=get(h_plot(kt),'xlabel');
    set(h_label,'fontsize',20);
    set(h_label,'fontweight','bold');
    h_label=get(h_plot(kt),'ylabel');
    set(h_label,'fontsize',20);
    set(h_label,'fontweight','bold');
    h_label=get(h_plot(kt),'zlabel');
    set(h_label,'fontsize',20);
    set(h_label,'fontweight','bold');    
    h_label=get(h_plot(kt),'title');
    set(h_label,'fontsize',20);
    set(h_label,'fontweight','bold');    
  end
  for it=1:length(h_objects),
    if strcmp(get(h_objects(it),'Type'),'line'),
      set(h_objects(it),'LineWidth',3);
       set(h_objects(it),'Markersize',16);
      set(h_objects(it),'Markersize',8);
    end
    if strcmp(get(h_objects(it),'Type'),'text')
      set(h_objects(it),'fontsize',16);
      set(h_objects(it),'fontweight','bold');
    end
  end 
  set(gca, 'Layer','top')
  grid on
end
