function [T_sup, m_sup, varout, t_sol] = ctrl_fun(t,T_air)
    persistent F;
    if isempty(F)
        F = evalin('base','P.CtrlFun');
    end
    [T_sup, m_sup, varout, t_sol] = feval(F,t,T_air);
end

