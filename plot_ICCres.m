function plot_ICCres(ICCres,cmax)

	if ~exist('cmax','var')
		cmax = 0.45;
	end

	figure
	imagesc(ICCres);
	axis image;
	colormap(jet(256));
	caxis([0 cmax]);
	colorbar;
	xlabel('residue index');
	ylabel('residue index');
	title('coarse-grained ICC map');
end