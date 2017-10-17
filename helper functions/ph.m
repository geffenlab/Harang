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
         end
         
         function imgsqz(dat, varargin)
             imagesc(squeeze(dat), varargin{:})
         end
         
         function pltsqz(dat, varargin)
             plot(squeeze(dat), varargin{:})
         end
         
     end

end