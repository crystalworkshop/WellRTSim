function [alpha_g, alpha_l, rho_g, rho_l, h_g, h_l, mu_l, T] = ...
            calculatePhaseProperties(p, h, state)
    % Inputs:
    %   p - Pressure [Pa]
    %   h - Mixture enthalpy [J/kg]
    %   state - Structure with interpolants in state.prop and state.deriv
    %
    % Outputs:
    %   alpha_g - Gas volume fraction [-]
    %   alpha_l - Liquid volume fraction [-]
    %   rho_g - Gas phase density [kg/m^3]
    %   rho_l - Liquid phase density [kg/m^3]
    %   h_g - Gas phase enthalpy [J/kg]
    %   h_l - Liquid phase enthalpy [J/kg]
    %   mu_l - Liquid phase viscosity [Pa·s]
    %   T - Temperature [K]

    origSize = size(p);
    pVec = p(:);
    hVec = h(:);
    n = numel(pVec);
    alpha_g = zeros(n, 1);
    alpha_l = ones(n, 1);
    rho_g   = zeros(n, 1);
    rho_l   = zeros(n, 1);
    h_g     = zeros(n, 1);
    h_l     = zeros(n, 1);
    mu_l    = zeros(n, 1);
    T       = zeros(n, 1);

    if ~state.interp
        P = pVec / 1e6; % MPa
        H = hVec / 1e3; % kJ/kg
        T = IAPWS_IF97('T_ph', P, H);
        x_v = IAPWS_IF97('x_ph', P, H);
        mu_l = IAPWS_IF97('mu_ph', P, H);

        xTransition = 1e-5;
        smoothWidth = max(xTransition / 5, eps);
        phaseWeight = 0.5 * (tanh((x_v - xTransition) ./ smoothWidth) + 1);
        phaseWeight = min(max(phaseWeight, 0), 1);

        maskTwo = x_v > 2e-16;
        maskLiq = ~maskTwo;

        if any(maskLiq)
            rho_l(maskLiq) = 1 ./ IAPWS_IF97('v_ph', P(maskLiq), H(maskLiq));
            h_l(maskLiq) = hVec(maskLiq);
        end

        if any(maskTwo)
            rho_g(maskTwo) = 1 ./ IAPWS_IF97('vV_p', P(maskTwo));
            h_g(maskTwo)   = IAPWS_IF97('hV_p', P(maskTwo)) * 1e3;
            h_l(maskTwo)   = IAPWS_IF97('hL_p', P(maskTwo)) * 1e3;
            rho_l(maskTwo) = 1 ./ IAPWS_IF97('vL_p', P(maskTwo));
            alpha_g(maskTwo) = x_v(maskTwo) ./ ...
                (x_v(maskTwo) + (1 - x_v(maskTwo)) .* (rho_g(maskTwo) ./ rho_l(maskTwo)));
            alpha_l(maskTwo) = 1 - alpha_g(maskTwo);
        end

        if any(maskTwo)
            alpha_g(maskTwo) = alpha_g(maskTwo) .* phaseWeight(maskTwo);
            alpha_l(maskTwo) = 1 - alpha_g(maskTwo);
        end

        alpha_l(maskLiq) = 1;
    else
        Pc = pVec;
        Hc = hVec;
        T = state.Temp(Hc, Pc);
        h_g = state.h_v(Pc);
        h_l = state.h_l(Pc);
        mu_l = state.visc(Hc, Pc);

        denom = h_g - h_l;
        x_v = (hVec - h_l) ./ denom;

        xTransition = 1e-5;
        smoothWidth = max(xTransition / 5, eps);
        phaseWeight = 0.5 * (tanh((x_v - xTransition) ./ smoothWidth) + 1);
        phaseWeight = min(max(phaseWeight, 0), 1);

        maskTwo = x_v > 1e-16;
        maskLiq = ~maskTwo;

        if any(maskLiq)
            rho_l(maskLiq) = state.dens(Hc(maskLiq), Pc(maskLiq));
            h_l(maskLiq) = hVec(maskLiq);
            h_g(maskLiq) = 0;
        end

        if any(maskTwo)
            rho_l(maskTwo) = state.dens_l(Pc(maskTwo));
            rho_g(maskTwo) = state.dens_v(Pc(maskTwo));
            alpha_g(maskTwo) = x_v(maskTwo) ./ ...
                (x_v(maskTwo) + (1 - x_v(maskTwo)) .* (rho_g(maskTwo) ./ rho_l(maskTwo)));
            alpha_l(maskTwo) = 1 - alpha_g(maskTwo);
        end

        if any(maskTwo)
            alpha_g(maskTwo) = alpha_g(maskTwo) .* phaseWeight(maskTwo);
            alpha_l(maskTwo) = 1 - alpha_g(maskTwo);
        end

        alpha_l(maskLiq) = 1;
    end

    alpha_g = reshape(alpha_g, origSize);
    alpha_l = reshape(alpha_l, origSize);
    rho_g   = reshape(rho_g,   origSize);
    rho_l   = reshape(rho_l,   origSize);
    h_g     = reshape(h_g,     origSize);
    h_l     = reshape(h_l,     origSize);
    mu_l    = reshape(mu_l,    origSize);
    T       = reshape(T,       origSize);
end
