classdef ph
    %PH Plot Helper
    %   has useful functions

     properties
     end

     methods (Static)
         
         function prefs()
             box off
             set(gca,'LineWidth',2);
             set(gca,'TickDir','out');
             set(gca,'FontSize',16);
             set(gca,'DefaultFigureWindowStyle','docked')
%              set(findall(gca, 'Type', 'Line'),'LineWidth',2);
         end
         
         function h = imgsqz(dat, varargin)
             h = imagesc(squeeze(dat), varargin{:});
         end
         
         function h = pltsqz(dat, varargin)
             h = plot(squeeze(dat), varargin{:});
         end
         
         function plot_ca_img(spatialInfo, activity)
             % activity tranlucency values, 0<=x<=1
             imagesc(spatialInfo.im)
             axis square
             for i = 1 : length(spatialInfo.ipix)
                 roi = spatialInfo.ROIs{i};
                 patch(roi(:,2),roi(:,1),'g','FaceAlpha',activity(i),...
                     'LineStyle','none');
             end
         end
         
         function error_shade(x, y, e, c, varargin)
             alpha = 0.2;
             plot(x,y,c,varargin{:})
             hold on
             patch([x flip(x)], [y-e flip(y+e)],c,'FaceColor',c,...
                 'FaceAlpha',alpha,'LineStyle','none')
             hold off
         end
         
     end

end