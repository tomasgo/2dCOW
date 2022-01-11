function [optim_pars,OS,diagnos] = optim_cow(X,optim_space,options,ref)
% function [optim_pars,OS,diagnos] = optim_cow(X,optim_space,options,ref)
% The routine automatically optimizes the segment length and slack size for COW
% alignment pre-processing.
% It does a "discrete-coordinates simplex" optimization (EVOP-like) for
% segment and slack parameters in COW alignment with so-called
% “Warping Effect = "Simplicity" + “Peak Factor”
% FvdB/ThS 061029
%
% in: X (n x m) data table with "n" objects and "m" variables
%     optim_space (1 x 4) optimization space - segment minimum, maximum and slack minimum, maximum
%     options (1 x 4) 1 : trigger plot and progress text
%                     2 : number of optimizations from grid maxima
%                     3 : maximum number of optimization steps
%                     4 : Fraction of maximal deviation from center in COW alignment
%                     default [0 3 50 0.15] (no plot; 3 starts from 3 maxima in grid search; maximum 50 steps; 15%)
%     ref (1 x m) reference object used in COW alignment (vector); if omitted
%                 reference is selected from the matrix "X" by "ref_select.m" with option 5
%
% out: optim_pars (1 x 2) optimal segment length and slack size
%      OS (7 x N) optimization sequence (first row segment, second
%      slack,third n. of added segments
%                                       third row "Warping Effect", fourth "Simplicity", Fifth "Peak Factor")
%      diagnos (struct): simplicity raw data, total run time, start points for optimization (columns in OS),
%                        "optim_space" and warping results for optimum (path + aligned matrix + diagnostics)
%
% uses ref_select.m, cow.m
%
% Authors:
% Thomas Skov / Frans van den Berg
% Royal Agricultural and Veterinary University - Department of Food Science
% Quality and Technology - Spectroscopy and Chemometrics group - Denmark
% email: thsk@kvl.dk / fb@kvl.dk - www.models.kvl.dk

if (nargin < 3)
   help optim_cow;
   return;
elseif (nargin == 2)
   options = [0 3 50 0.15];
elseif (nargin == 3)
   [ref,~,refN] = ref_select(X,[],[5 0]);
   if options(1)
      disp(['Object ' num2str(refN) ' selected as reference']);
   end
end

if (length(options) ~= 4)
   options = [0 3 50 0.15];
end

[Nspace,Mspace] = size(optim_space);
if (Nspace ~= 1) || (Mspace ~= 4)
   error('ERROR: "optim_space" must be of size (1 x 4)');
end

S  = sum(svd(X/sqrt(sum(X(:).^2))).^4);
ag = unique(linspace(optim_space(1),optim_space(2)),4);
bg = unique(linspace(optim_space(3),optim_space(4),4));

t00 = clock;
if options(1)
   disp('Starting grid search');
end
N        = 1;
nGridRun = length(ag) * length(bg);
OS       = zeros(5,nGridRun);
OS(1,:)  = a(ones(length(bg),1),:);
OS(2,:)  = b(ones(length(ag),1),:)';

for i_grid = 1:nGridRun

   if options(1)
      t0 = clock;
   end
   [OS(:,i_grid),exitflag] = optim_eval(X,OS(1:2,i_grid),OS,ref,round(length(ref)*options(4))); pause(5);
   if options(1)
      if (exitflag==1)
         s = fprintf('\nrun %i/%i: - min (segment/slack combination was already computed)',i_grid,nGridRun);
      elseif (exitflag==2)
         s = fprintf('\nrun %i/%i: - min (illegal segment/slack combination)',i_grid,nGridRun);
      else
         s = fprintf('\nrun %i/%i: %3.2f min',i_grid,nGridRun,etime(clock,t0)/60);
      end
      
   end

end

[~,c] = unique(OS(3,:));
starts    = fliplr(c(end - options(2) + 1:end));
N         = nGridRun + 1;
for a = 1:length(starts)

   if options(1)
      disp(['Starting optimization ' num2str(a) '/' num2str(length(starts)) ', (segment = ' num2str(OS(1,starts(a))) ', slack = ' num2str(OS(2,starts(a))) ')']);
      t0 = clock;
   end
   Na        = N - 1;
   ps        = [starts(a) 0 0];
   
   OS(1:2,N) = OS(1:2,ps(1)) + [1 0]';
   [OS(:,N)] = optim_eval(X,OS(1:2,N),OS,ref,round(length(ref)*options(4))); pause(5);

   ps(2)     = N;
   N         = N + 1;
   OS(1:2,N) = OS(1:2,ps(1)) + [0 1]';
   [OS(:,N)] = optim_eval(X,OS(1:2,N),OS,ref,round(length(ref)*options(4))); pause(5);
   ps(3)     = N;
   N         = N + 1;

   pt = 1;
   while pt
      
      [b,c]     = sort(OS(3,ps));
      OS(1:2,N) = OS(1:2,ps(c(1))) + sum(OS(1:2,ps)>OS(1:2,ps(c(1))),2) - sum(OS(1:2,ps)<OS(1:2,ps(c(1))),2);
      [OS(:,N)] = optim_eval(X,OS(1:2,N),OS,ref,round(length(ref)*options(4))); pause(5);
      if (OS(3,N) <= b(1))

         N         = N + 1;
         c         = c([2 1 3]);
         b         = b([2 1 3]);
         OS(1:2,N) = OS(1:2,ps(c(1))) + sum(OS(1:2,ps)>OS(1:2,ps(c(1))),2) - sum(OS(1:2,ps)<OS(1:2,ps(c(1))),2);
         [OS(:,N)] = optim_eval(X,OS(1:2,N),OS,ref,round(length(ref)*options(4))); pause(5);
         if (OS(3,N) <= b(1))
            pt = 0;
         else
            ps(c(1)) = N;
         end

      end
      ps(c(1)) = N;
      N = N + 1;
      if ((N - Na - 1) >= options(3))
         pt = 0;
         disp(['   Early termination after ' num2str(N-Na-1) ' steps!']);
      end

   end
   if options(1)
      s = ['optimization ' num2str(a) '/' num2str(length(starts)) ' : ' num2str(etime(clock,t0)/60,2) 'min, ' num2str(N-Na-1) ' steps'];
      disp(s);
   end
   steps(a) = N - Na-1;

end

[~,optim] = max(OS(3,:));
optim_pars    = OS(1:2,optim);

% Plotting
if options(1)
   f = figure;
   stem3(OS(1,:),OS(2,:),OS(3,:),'filled');
   xlabel('Segment length'); ylabel('Slack size'); zlabel('Warping Effect');
   s = ['Warping Effect(' num2str(optim_pars(1)) ',' num2str(optim_pars(2)) ')'];
   text(OS(1,optim),OS(2,optim),OS(3,optim),s);

   figure;
   stem3(OS(1,:),OS(2,:),OS(4,:),'filled');
   hold on;
   stem3(OS(1,:),OS(2,:),ones(size((OS(4,:))))*S);
   hold off;
   title('(Solid = Simplicity, Open = Simplicity Raw Data)');
   xlabel('Segment length'); ylabel('Slack size'); zlabel('Simplicity');
   [~,b] = max(OS(4,:));
   s = ['Simplicity(' num2str(OS(1,b)) ',' num2str(OS(2,b)) ')'];
   text(OS(1,b),OS(2,b),OS(4,b),s);
   s = ['Warping Effect(' num2str(optim_pars(1)) ',' num2str(optim_pars(2)) ')'];
   text(OS(1,optim),OS(2,optim),OS(4,optim),s);

   figure;
   stem3(OS(1,:),OS(2,:),OS(5,:),'filled');
   xlabel('Segment length'); ylabel('Slack size'); zlabel('Peak Factor');
   [~,b] = max(OS(5,:));
   s = ['Peak Factor(' num2str(OS(1,b)) ',' num2str(OS(2,b)) ')'];
   text(OS(1,b),OS(2,b),OS(5,b),s);
   s = ['Warping Effect(' num2str(optim_pars(1)) ',' num2str(optim_pars(2)) ')'];
   text(OS(1,optim),OS(2,optim),OS(5,optim),s);
   figure(f);
   disp(['Finished optimization, optimal (segment slack) = (' num2str(optim_pars') '), total time : ' num2str(etime(clock,t00)/60,2) 'min']);
end

% Diagnostics
if (nargout > 2)
   if options(1)
      disp('Computing diagnostics');
   end
   diagnos.base_simplicity = S;
   diagnos.time_min = etime(clock,t00)/60;
   diagnos.optim_starts_in_OS = starts;
   diagnos.optim_steps_in_OS = steps;
   diagnos.optim_space = optim_space;
   diagnos.reference = ref;
   if exist('refN','var')
      diagnos.reference_sample = refN;
   end
   try
      [diagnos.warping,diagnos.Xw,diagnos.warping_diagnos] = cow(ref,X,optim_pars(1),optim_pars(2),[0 1 0 round(length(ref)*options(4)) 0]);
   catch
      disp('Data is to big for "fast COW" implementation, final result not included in "diagnos".');
   end
   if options(1)
      figure
      subplot(2,1,1); plot(1:size(X,2),X);
      title('Data raw');
      grid;
      subplot(2,1,2); plot(1:size(diagnos.Xw,2),diagnos.Xw);
      title(['Data from optimal correction (segment ' num2str(optim_pars(1)) ', slack ' num2str(optim_pars(2)) ')']);
      grid;
   end
end

%%%
function [y,exitflag] = optim_eval(X,p,OS,ref,losange)
index1 = find(OS(1,1:end-1)==p(1));
index2 = find(OS(2,1:end-1)==p(2));
index3 = intersect(index1,index2);
y(1:2) = p;
exitflag = 0;
if index3
   y(3) = OS(3,index3(1));
   y(4) = OS(4,index3(1));
   y(5) = OS(5,index3(1));
   exitflag = 1;
else
   if (p(1) <= p(2)+3) || (p(2) < 1)% segment > slack OR slack < 1
      y(3:5) = 0;
      exitflag = 2;
   else
      K = size(X,1);
      for k = 1:K
         normX(k) = norm(X(k,:));
      end
      try
         [warping,X,diagnos] = cow(ref,X,p(1),p(2),[0 1 0 losange 0]);
      catch
         for a=1:K
            [warping,X(a,:),diagnos] = cow(ref,X(a,:),p(1),p(2),[0 1 0 losange 0]);
         end
      end
      y(1) = diagnos.segment_length(1,1)+1;
      y(2) = diagnos.slack;
      y(4) = (sum(svds(X/sqrt(sum(X(:).^2)),K).^4));
      for k = 1:K
         PEAKFAC(k) = abs(((norm(X(k,:))-normX(k)))/normX(k));
         PEAKFAC(k) = (1-min([PEAKFAC(k) 1]))^2;
      end
      y(5) = mean(PEAKFAC);
      y(3) = y(4) + y(5);
   end
end