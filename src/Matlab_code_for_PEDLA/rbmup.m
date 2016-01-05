function x = rbmup(rbm, x)
    x = sigm(bsxfun(@plus, rbm.c', x * rbm.W'));
end
