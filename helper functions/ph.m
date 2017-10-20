classdef ph
    %PH Plot Helper
    %   has useful functions

     properties
     end

     methods (Static)
         
         function prefs
             box off
             set(gca,'LineWidth',2);
             set(gca,'TickDir','out');
             set(gca,'FontSize',16);
             set(gca,'DefaultFigureWindowStyle','docked')
             set(findall(gca, 'Type', 'Line'),'LineWidth',2);
         end
         
         function h = imgsqz(dat, varargin)
             h = imagesc(squeeze(dat), varargin{:});
         end
         
         function h = pltsqz(dat, varargin)
             h = plot(squeeze(dat), varargin{:});
         end
         
     end

end