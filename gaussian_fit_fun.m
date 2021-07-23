%{
Author: S. Watanabe
Last updated; 7/22/2021

For least square Gaussian fitting by the model
  data = offset + amp * exp(-(x - mu)^2 / (2 * sigma^2))

Input
  data - vector of data to fit
  parameters - initial values of fitting (optional)

Output
  params - structure of offset, amp, mu, sigma
  sse - mean square error of fitting

Usage
  [params, sse] = gaussian_fit_fun(ydata, xdata, 'offset', offset_init, 'amp', amp_init, 'mu', mu_init, 'sigma', sigma_init)
  [params, sse] = gaussian_fit_fun(ydata, 'offset', offset_init, 'amp', amp_init, 'mu', mu_init, 'sigma', sigma_init)

  When xdata is omitted, xdata is assigned as 1:length(ydata).
  When some parameters are omitted, default values are used (see the program below).

%}


function [params, sse] = gaussian_fit_fun(data, varargin)
	nvarargin = length(varargin);
	
	nlinept = length(data);

	xdata = (1:nlinept)';
	ydata = data(:);
	ydatafit = zeros(nlinept, 1);

	%default initial values
	offset_init = min(ydata);
	[amp_init, mu_init] = max(ydata - offset_init);
	sigma_init = sum(ydata > offset_init + amp_init/2) /2;

	for i=1:nvarargin
		if isnumeric(varargin{i}) && i==1
			xdata = varargin{i}(:);
		end
		
		if ischar(varargin{i})
			switch varargin{i}
				case 'offset'
					offset_init = varargin{i+1};
				case 'amp'
					amp_init = varargin{i+1};
				case 'mu'
					mu_init = varargin{i+1};
				case 'sigma'
					sigma_init = varargin{i+1};
			end
		end
	end
	
	params_init = [offset_init, amp_init, mu_init, sigma_init];
	[params, sse] = fminsearch(@gaussian_fit_f, params_init);


	function sse = gaussian_fit_f(params)
		offset = params(1);
		amp = params(2);
		mu = params(3);
		sigma = params(4);

		ydatafit = offset + amp * exp(-(xdata - mu).^2 / (2 * sigma^2) );
		sse = mean( (ydatafit - ydata).^2 );
	end

end
