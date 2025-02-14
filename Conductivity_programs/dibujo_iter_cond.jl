using JLD, Plots, Measures, LaTeXStrings

err_l2_s = FileIO.load("results_iter_s.jld2","err_l2")
err_lin_s = FileIO.load("results_iter_s.jld2","err_lin")
err_l2_d = FileIO.load("results_iter_d.jld2","err_l2")
err_lin_d = FileIO.load("results_iter_d.jld2","err_lin")
iter1=0 : length(err_l2_s)-1
iter2=0 : length(err_l2_d)-1

#p1=scatter(iter2,log10.(err_l2_s))   
#p1=scatter!(iter2,log10.(err_l2_d)) 
p1= plot(iter2,log10.(err_l2_s),label = L"q_{smooth}",linewidth=2,marker=(:circle,5),
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     xlabel = "iterations",ylabel = "\$ log_{10}(||q-q^n||_{L^2(0,1)}) \$", bottom_margin=5mm, right_margin=10mm, left_margin=10mm)
  
p1= plot!(iter2,log10.(err_l2_d),label = L"q_{disc}",legendfontsize=16,linewidth=2,marker=(:square,5))

display(p1)
savefig(p1,"ejem_ga_iter_conv_l2.png")

p2= plot(iter2,log10.(err_lin_s),label = L"q_{smooth}",linewidth=2,marker=(:circle,5),ylabel = L"log_{10}(||q-q^n||_{L^\infty(0,1)})",
     xtickfontsize=16,ytickfontsize=16,legendfontsize=16,xguidefontsize=16,yguidefontsize=16,
     xlabel = "iterations", bottom_margin=5mm, right_margin=10mm, left_margin=10mm)
  
p2= plot!(iter2,log10.(err_lin_d),label = L"q_{disc}",legendfontsize=16,linewidth=2,marker=(:square),markersize=5) # ,ylims=(-2.1,0.5))

display(p2)
savefig(p2,"ejem_ga_iter_conv_lin.png")