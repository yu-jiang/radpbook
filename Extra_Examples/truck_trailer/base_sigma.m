%%
function s = base_sigma(x1, x2, x3)
    s = [x1 * (x2 - sin(x2))
         x2 * (x2 - sin(x2))
         x3 * (x2 - sin(x2))
         x1 * (x3 - sin(x3))
         x2 * (x3 - sin(x3))
         x3 * (x3 - sin(x3))
         ];       
end