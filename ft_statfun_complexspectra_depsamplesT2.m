function [s, cfg] = ft_statfun_complexspectra_depsamplesT2(cfg, dat, design)

% FT_STATFUN_COMPLEXSPECTRA_DEPSAMPLES calculates the dependent samples Hotelling's T-squared test  
% on the biological complex spectra data in dat (the dependent real and imaginary variables), 
% using the information on the independent variable (ivar) in design.
%
% Use this function by calling the high-level statistics functions 
%   [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
% with the following configuration option
%   cfg.statistic = 'complexspectra_depsamplesT2'
% see FT_FREQSTATISTICS for details.
%
% For low-level use, the external interface of this function has to be
%   [s,cfg] = statfun_complexspectra_depsamplesT2(cfg, dat, design);
% where
%   dat    contains the biological complex spectra data, Nsamples
%   design contains the independent variable (ivar) and the unit-of-observation (uvar) 
%          factor,  Nfac
%
% Configuration options
%   cfg.computestat    = 'yes' or 'no', calculate the statistic (default='yes')
%   cfg.computecritval = 'yes' or 'no', calculate the critical values of the test statistics (default='no')
%   cfg.computeprob    = 'yes' or 'no', calculate the p-values (default='no')
% The following options are relevant if cfg.computecritval='yes' and/or
% cfg.computeprob='yes'.
%   cfg.alpha = critical alpha-level of the statistical test (default=0.05)
%   cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%               cfg.tail in combination with cfg.computecritval='yes'
%               determines whether the critical value is computed at
%               quantile cfg.alpha (with cfg.tail=-1), at quantiles
%               cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%               quantile (1-cfg.alpha) (with cfg.tail=1).
%
% Design specification
%   cfg.ivar  = row number of the design that contains the labels of the conditions that must be 
%               compared (default=1). The labels are the numbers 1 and 2.
%   cfg.uvar  = row number of design that contains the labels of the units-of-observation (subjects or trials)
%               (default=2). The labels are assumed to be integers ranging from 1 to 
%               the number of units-of-observation.
% 
% (The T2 calculations are implemented in modified versions of FT_STATFUN_DEPSAMPLEST and FT_STATFUN_DEPSAMPLESF 
% by Niels Trusbak Haumann, Center for Music in the Brain, Aarhus University, Denmark, February 2024, 
% with inspiration from A. Trujillo-Ortiz & R. Hernandez-Walls, Copyright (C) December 2002 and Picton et al. (2003).) 
% 
% This file is compatible with FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

% set defaults
if ~isfield(cfg, 'computestat'),    cfg.computestat    = 'yes'; end 
if ~isfield(cfg, 'computecritval'), cfg.computecritval = 'no';  end
if ~isfield(cfg, 'computeprob'),    cfg.computeprob    = 'no';  end
if ~isfield(cfg, 'alpha'),          cfg.alpha          = 0.05;  end
if ~isfield(cfg, 'tail'),           cfg.tail           = 1;     end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end;
if ~isfield(cfg,'uvar') || isempty(cfg.uvar)
    error('uvar must be specified for dependent samples statistics');
end

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
n1  = length(sel1);
n2  = length(sel2);
if (n1+n2)<size(design,2) || (n1~=n2)
  error('Invalid specification of the design array.');
end
nunits = length(design(cfg.uvar, sel1));
if nunits<2
    error('The data must contain at least two units (usually subjects).')
end
if (nunits*2)~=(n1+n2)
  error('Invalid specification of the design array.');
end
nsmpls = size(dat,1);

if strcmp(cfg.computestat,'yes')
  % compute the statistic
  % store the positions of the 1-labels and the 2-labels in a nunits-by-2 array
  poslabelsperunit = zeros(nunits,2);
  poslabel1        = find(design(cfg.ivar,:)==1);
  poslabel2        = find(design(cfg.ivar,:)==2);
  [dum,i]          = sort(design(cfg.uvar,poslabel1), 'ascend');
  poslabelsperunit(:,1) = poslabel1(i);
  [dum,i]          = sort(design(cfg.uvar,poslabel2), 'ascend');
  poslabelsperunit(:,2) = poslabel2(i);
    
  % calculate the differences between the conditions
  diffmat = zeros(nsmpls,nunits,2); % Difference matrix ( channel-freq-time vector x subjects x real/imag axes )
  diffmat(:,:,1) = real(dat(:,poslabelsperunit(:,1))) - real(dat(:,poslabelsperunit(:,2)));
  diffmat(:,:,2) = imag(dat(:,poslabelsperunit(:,1))) - imag(dat(:,poslabelsperunit(:,2)));
    
  % Calculate the dependent samples Hotelling's T-Squared statistics
  avgdiff = nanmean(diffmat,2); % Average difference vectors
  diffmatc = diffmat - repmat(avgdiff,[1,size(diffmat,2),1]); % Centered difference matrix
  covmatdiff = []; % Variance-covariance matrix ( channel-freq-time vector x real/imag axes )...
  covmatdiff(:,1,1) = 1/(nunits-1) * nansum( diffmatc(:,:,1).^2 , 2 );
  covmatdiff(:,2,2) = 1/(nunits-1) * nansum( diffmatc(:,:,2).^2 , 2 );
  covmatdiff(:,1,2) = 1/(nunits-1) * nansum( diffmatc(:,:,1).*diffmatc(:,:,2) , 2 );
  covmatdiff(:,2,1) = covmatdiff(:,1,2);
  inv_det = 1./( covmatdiff(:,1,1).*covmatdiff(:,2,2) - covmatdiff(:,1,2).*covmatdiff(:,2,1)); % Inverse determinant for variance-covariance matrix
  covmatdiff_inv = zeros(size(covmatdiff)); % Inverse variance-covariance matrix...
  invid = find(1./inv_det~=0); % Use the fast standard inverse calculation when the determinant differs from zero
  if ~isempty(invid)
      covmatdiff_inv(invid,1,1) = inv_det(invid) .* covmatdiff(invid,2,2);
      covmatdiff_inv(invid,1,2) = inv_det(invid) .* -covmatdiff(invid,1,2);
      covmatdiff_inv(invid,2,1) = inv_det(invid) .* -covmatdiff(invid,2,1);
      covmatdiff_inv(invid,2,2) = inv_det(invid) .* covmatdiff(invid,1,1);
  end
  pinvid = find(1./inv_det==0); % Use the More-Penrose pseudo-inverse estimate when the determinant is zero
  if ~isempty(pinvid)
      for i=1:length(pinvid)
          covmatdiff_inv(pinvid(i),:,:) = pinv(squeeze(covmatdiff(pinvid(i),:,:)));
      end
  end
  s.stat = nunits * sum( squeeze( sum( repmat( squeeze(avgdiff), [1, 1, 2]) .* covmatdiff_inv, 2) ) .* squeeze(avgdiff), 2); % Hotelling's T-Squared statistic
end

if strcmp(cfg.computecritval,'yes')
  % Also compute the critical values
  s.v1      = 2;        % Numerator degrees of freedom
  s.v2      = nunits-2; % Denominator degrees of freedom
  if cfg.tail==-1
      error('For a dependent samples T-squared-statistic, it does not make sense to calculate a left tail critical value.');
  end;
  if cfg.tail==0
      error('For a dependent samples T-squared-statistic, it does not make sense to calculate a two-sided critical value.');
  end;
  if cfg.tail==1
      if nunits < 50 % F-distribution
          s.critval = finv(1-cfg.alpha, s.v1, s.v2);
      elseif nunits >= 50 % Chi-square approximation
          s.critval = chi2inv(1-cfg.alpha, 2);
      end
  end
end

if strcmp(cfg.computeprob,'yes')
  % Also compute the p-values
  s.v1      = 2;        % Numerator degrees of freedom
  s.v2      = nunits-2; % Denominator degrees of freedom
  F = (nunits-2)/(2*nunits-2)*s.stat;
  if cfg.tail==-1
      error('For a dependent samples T-squared-statistic, it does not make sense to calculate a left tail critical value.');
  end;
  if cfg.tail==0
      error('For a dependent samples T-squared-statistic, it does not make sense to calculate a two-sided critical value.');
  end;
  if cfg.tail==1
      if nunits < 50 % F-distribution
          s.prob = 1-fcdf(F, s.v1, s.v2);
      elseif nunits >= 50 % Chi-square approximation
          s.prob = 1-chi2cdf(s.stat, 2);
      end
  end
end

