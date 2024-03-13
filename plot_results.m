(* ::Package:: *)

(* ::Input:: *)
(*data= Import[File["C:\\c++ proj\\flow_dynamics\\flow_dynamic\\build\\results\\log_integrals.txt"], "Table"];*)
(*time = 100;*)
(*steps = 50;*)
(*tau = time / steps;*)
(*nu = 0.001;*)


(* ::Input:: *)
(*(*ResourceFunction["MaTeXInstall"][]*)
(*ConfigureMaTeX[*)
(*"pdfLaTeX"\[Rule]"C:\\texlive\\2023\\bin\\windows\\xelatex.exe","Ghostscript"\[Rule]"C:\\Program Files\\gs\\gs10.03.0\\bin\\gswin64c.exe"*)
(*]*)
(*<<MaTeX`*)*)
(**)


(* ::Input:: *)
(*Grid[data,Frame->All];*)
(*omega = data[[3;;,1]];*)
(*omegasq= data[[3;;,2]];*)
(*epsilon = data[[3;;,3]];*)
(*domegadtpsi = data[[3;;,4]];*)
(*error = data[[3;;,5]];*)
(*depsilondt = Table[(epsilon[[i]] - epsilon[[i-1]])/tau,{i,2,Length@epsilon}];*)
(*domegadtpsi2 =  Table[(domegadtpsi[[i]] + domegadtpsi[[i-1]])/2,{i,2,Length@domegadtpsi}];*)
(*omegasqnu= Table[nu*(omegasq[[i]] + omegasq[[i-1]])/2,{i,2,Length@omegasq}];*)


(* ::Input:: *)
(*(*ListPlot[{omega},PlotLegends\[Rule]{"\[Omega]"}]*)
(*ListPlot[{omegasq},PlotLegends\[Rule]{"\[Omega]^2"}]*)
(*ListPlot[{epsilon},PlotLegends\[Rule]{"\[Epsilon]"}]*)
(*ListPlot[{domegadtpsi},PlotLegends\[Rule]{"d\[Omega]/dt\[Psi]"}]*)*)
(*ListPlot[{depsilondt,domegadtpsi2}, PlotLegends->{"\!\(\*FractionBox[\(d\[Epsilon]\), \(dt\)]\)","\[Integral]\!\(\*FractionBox[\(dw\), \(dt\)]\)\[Psi]\[DifferentialD]V"},PlotRange->Full]*)
(*ListLogPlot[{Abs[depsilondt],Abs[domegadtpsi2]}, PlotLegends->{"\!\(\*FractionBox[\(d\[Epsilon]\), \(dt\)]\)","\[Integral]\!\(\*FractionBox[\(dw\), \(dt\)]\)\[Psi]\[DifferentialD]V"},PlotRange->Full]*)
(*ListLogPlot[{Abs[#]&/@(depsilondt-domegadtpsi2)}, PlotLegends->{"\!\(\*FractionBox[\(d\[Epsilon]\), \(dt\)]\) -\[Integral] \!\(\*FractionBox[\(dw\), \(dt\)]\)\[Psi]\[DifferentialD]V"},PlotRange->Full]*)
(*(*ListPlot[{depsilondt,omegasqnu}, PlotLegends\[Rule]{"d\[Epsilon]/dt","\[Nu]\[Integral]\[Omega]^2\[DifferentialD]V"}]*)
(*ListLogPlot[{Abs[#]&/@(depsilondt-omegasqnu)}, PlotLegends\[Rule]{"d\[Epsilon]/dt -\[Nu]\[Integral] \[Omega]^2\[DifferentialD]V"}]*)*)
(*ListLogPlot[error]*)



