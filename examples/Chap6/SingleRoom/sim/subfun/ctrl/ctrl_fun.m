function [H_flo, H_win, varout] = ctrl_fun(t,T_air)
    persistent F;
    if isempty(F)
        F = evalin('base','P.CtrlFun');
    end
    [H_flo, H_win, varout] = feval(F,t,T_air);
end