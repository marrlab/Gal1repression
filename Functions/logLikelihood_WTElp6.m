function varargout = logLikelihood_WTElp6(xi,varargin)

n_xi = numel(xi);
options.ami.sensi_meth = 'forward';
options.ami.atol = 1e-15;
options.ami.rtol = 1e-8;
options.ami.sensi = 0;
%options.llh.scale = 'lin'; % data on log or linear scale, mainly needed for
%plotting reasons
% options.llh.distribution = 'normal'; % 'laplace'

DA = varargin{1};
imodel = varargin{2};

if imodel == 1
    indA = 1:4;
else
    indA = 1:5;
end

D.Y = DA(1).y; %artificial data matrix // needed for windows use

try
    nderiv = nargout-1;
    logL = 0;
    if(isfield(options,'sens_ind'))
        if(nderiv>=1)
            dlogL = zeros(length(options.sens_ind),1);
            options.ami.sensi = 1;
        end
    else
        if(nderiv>=1)
            dlogL = zeros(n_xi,1);
            options.ami.sensi = 1;
        end
    end
    
    %% Simulation
    
    par = 10.^(xi);

    if imodel == 1
        P01 = par(indA(1));
        b1 = par(indA(2));
        c1 = par(indA(3));
        sigmayA = par(indA(4))*ones(length(DA.y),1);
    else
        P01 = par(indA(1));
        t_rep1 = par(indA(2));
        b1 = par(indA(3));
        c1 = par(indA(4));
        sigmayA = par(indA(5))*ones(length(DA.y),1);
    end
    
    %WT simulation
    count = 1;
    for t = DA(1).t
       if imodel == 1
            f1(count) = b1/(c1)+(P01-b1/(c1))*exp(-c1*t);
        else
            if t<t_rep1
                f1(count) = b1/(c1)+(P01-b1/(c1))*exp(-c1*t);
            else
                P0_init = b1/(c1)+(P01-b1/(c1))*exp(-c1*t_rep1);
                f1(count) = P0_init*exp(-c1*(t-t_rep1));
            end
        end
        count = count+1;
    end
    f1 = f1';
    
    %% Likelihood evaluation
    
    for idata = 1:length(DA)
        logL = logL - 0.5*sum(sum(~isnan(DA(idata).y)))*log(2*pi)-...
            nansum(nansum(log(sigmayA(~isnan(DA(idata).y))))) -...
            0.5*nansum(nansum(((f1-DA(idata).y)./sigmayA).^2));
    end
    
    %% Output assignment
    varargout{1} = logL;
    
catch
    varargout{1} = NaN;
end


