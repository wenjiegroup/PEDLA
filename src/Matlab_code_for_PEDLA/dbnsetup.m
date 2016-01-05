function dbn = dbnsetup(sizes, x, opts)
    n = size(x, 2);
	
	%确保网络结构为一个行向量
	if size(sizes,1)~=1 && size(sizes,2)~=1 
		error('sizes must be a vector!');
	elseif size(sizes,2)==1
		sizes=sizes';
	end
    dbn.sizes = [n, sizes];

    for u = 1 : numel(dbn.sizes) - 1
        dbn.rbm{u}.alpha    = opts.alpha;
        dbn.rbm{u}.momentum = opts.momentum;

        dbn.rbm{u}.W  = zeros(dbn.sizes(u + 1), dbn.sizes(u));
        dbn.rbm{u}.vW = zeros(dbn.sizes(u + 1), dbn.sizes(u));

        dbn.rbm{u}.b  = zeros(dbn.sizes(u), 1);
        dbn.rbm{u}.vb = zeros(dbn.sizes(u), 1);

        dbn.rbm{u}.c  = zeros(dbn.sizes(u + 1), 1);
        dbn.rbm{u}.vc = zeros(dbn.sizes(u + 1), 1);
    end

end
