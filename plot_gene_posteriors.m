function [fig, t, axs] = plot_gene_posteriors(gene, ind, ABC, SFABC, ks_true, varargin)
%PLOT_GENE_POSTERIORS Plot the 3-by-3 posterior block for one gene.
%
%   [fig,t,axs] = plot_gene_posteriors(gene,ind,ABC,SFABC,ks_true,...)
%
%   Inputs:
%     gene    - Gene index.
%     ind     - SMC iteration index.
%     ABC     - Structure with ks_t and ws_t; use [] to plot SFABC alone.
%     SFABC   - Structure with the same fields.
%     ks_true - Number-of-genes-by-3 array of true parameter values;
%               use [] to omit the true-parameter markers.
%
%   Array shapes (parameters must be on the original, positive scale):
%     ks_t: [number of genes, number of iterations, number of particles, 3]
%     ws_t: [number of genes, number of iterations, number of particles]
%     ABC and SFABC may contain different numbers of particles.
%
%   Name-value options:
%     'N'             - Resampled particles per method (default 1e5).
%     'NumBins'       - Bins for marginal histograms (default 25).
%     'NumBoxes'      - Bins per axis for joint histograms (default 25).
%     'Limits'        - Histogram limits in log10 space (default [-1 1]).
%     'TrueColor'     - RGB colour of truth markers (default orange).
%     'ABCColormap'   - Colormap used as supplied (default flipud(pink(256))).
%     'SFABCColormap' - Colormap used as supplied (default flipud(bone(256))).
%     'Layout'        - Existing tiledlayout; default creates a new figure.
%     'Tiles'         - 3-by-3 tile index matrix within Layout.
%                       Default: [1 2 3; 4 5 6; 7 8 9].
%
%   Outputs:
%     fig      - Figure handle; export with exportgraphics.
%     t        - Tiled layout handle.
%     axs(r,c) - Axes handle at block row r, column c.
%
%   Notes:
%     Diagonal: both marginals. Lower triangle: ABC. Upper triangle: SFABC.
%     With ABC = [], SFABC occupies the diagonal and upper triangle.
%     Call rng(...) before plotting for reproducible resampling.
%     Requires randsample (Statistics and Machine Learning Toolbox).
%
%   Example (run from the results folder):
%     ABC = load('ABC_SMC.mat');
%     SFABC = load('SFABC_SMC.mat');
%     D = load('Data.mat','ks');
%     rng(314)
%     fig = plot_gene_posteriors(1,numel(SFABC.tols),ABC,SFABC,D.ks);
%
%   Example after running SFABC inference in the workspace:
%     SFABC = struct('ks_t',ks_t,'ws_t',ws_t);
%     fig = plot_gene_posteriors(1,numel(tols),[],SFABC,ks);

p = inputParser;
isPosInt = @(x) isnumeric(x) && isscalar(x) && isreal(x) && ...
    isfinite(x) && x >= 1 && x == fix(x);
isMap = @(x) isnumeric(x) && isreal(x) && ismatrix(x) && ...
    ~isempty(x) && size(x,2) == 3 && ...
    all(isfinite(x(:))) && all(x(:) >= 0 & x(:) <= 1);
addParameter(p, 'N', 1e5, isPosInt);
addParameter(p, 'NumBins', 25, isPosInt);
addParameter(p, 'NumBoxes', 25, isPosInt);
addParameter(p, 'Limits', [-1 1], @(x) isnumeric(x) && isreal(x) && ...
    numel(x) == 2 && all(isfinite(x(:))) && x(1) < x(2));
addParameter(p, 'TrueColor', [255 164 0]/255, ...
    @(x) isMap(x) && size(x,1) == 1);
addParameter(p, 'ABCColormap', flipud(pink(256)), isMap);
addParameter(p, 'SFABCColormap', flipud(bone(256)), isMap);
addParameter(p, 'Layout', []);
addParameter(p, 'Tiles', [1 2 3; 4 5 6; 7 8 9], ...
    @(x) isnumeric(x) && isreal(x) && isequal(size(x),[3 3]) && ...
    all(isfinite(x(:))) && all(x(:) >= 1 & x(:) == fix(x(:))) && ...
    numel(unique(x(:))) == 9);
parse(p, varargin{:});
o = p.Results;
lims = reshape(o.Limits, 1, 2);

validateattributes(gene, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(ind, {'numeric'}, {'scalar','integer','positive','finite'});

% Select each gene/iteration, resample by its weights, then take log10.
hasABC = ~isempty(ABC);
ABC_ks = [];
if hasABC
    ABC_ks = resamplePosterior(ABC, gene, ind, o.N);
end
SFABC_ks = resamplePosterior(SFABC, gene, ind, o.N);
truth = [];
if ~isempty(ks_true)
    if size(ks_true,2) ~= 3 || gene > size(ks_true,1)
        error('ks_true must have three columns and include the selected gene.');
    end
    validateattributes(ks_true(gene,:), {'numeric'}, ...
        {'real','finite','positive'});
    truth = log10(ks_true(gene,:));
end

if isempty(o.Layout)
    fig = figure('Units','centimeters','Position',[5 2 12 12]);
    t = tiledlayout(fig, 3, 3, 'TileSpacing','tight');
else
    t = o.Layout;
    fig = ancestor(t, 'figure');
end
axs = gobjects(3,3);

for row = 1:3
    for col = 1:3
        ax = nexttile(t, o.Tiles(row,col));
        axs(row,col) = ax;
        if ~hasABC && row > col
            axis(ax, 'off');
            continue
        end
        hold(ax, 'on');

        if row == col
            % Marginals of both methods share the same bins.
            if hasABC
                histogram(ax, ABC_ks(row,:), 'LineWidth',1.25, ...
                    'EdgeColor','#B49D9B', 'LineStyle','-', ...
                    'DisplayStyle','stairs', 'Normalization','probability', ...
                    'NumBins',o.NumBins, 'BinLimits',lims);
            end
            histogram(ax, SFABC_ks(row,:), 'LineWidth',1.25, ...
                'EdgeColor','#A5C4C5', 'LineStyle','-', ...
                'DisplayStyle','stairs', 'Normalization','probability', ...
                'NumBins',o.NumBins, 'BinLimits',lims);
            xlabel(ax, sprintf('log_{10} k_%d', row));
            ylabel(ax, 'Frequency');
            xlim(ax, lims);
        else
            % Both triangles use x = lower-index parameter, y = higher.
            par1 = min(row,col);
            par2 = max(row,col);
            if hasABC && row > col
                samples = ABC_ks;
                cmap = o.ABCColormap;
                method = 'ABC';
            else
                samples = SFABC_ks;
                cmap = o.SFABCColormap;
                method = 'SFABC';
            end
            histogram2(ax, samples(par1,:), samples(par2,:), ...
                [o.NumBoxes o.NumBoxes], 'XBinLimits',lims, ...
                'YBinLimits',lims, 'Normalization','probability', ...
                'DisplayStyle','tile', 'ShowEmptyBins','on', 'EdgeAlpha',0);
            colormap(ax, cmap);
            title(ax, method, 'FontName','Arial', 'FontSize',7);
            xlabel(ax, sprintf('log_{10} k_%d', par1));
            ylabel(ax, sprintf('log_{10} k_%d', par2));
            xlim(ax, lims + [-0.1 0.1]);
            ylim(ax, lims + [-0.1 0.1]);
            ax.SortMethod = 'childorder';
            if ~isempty(truth)
                hTrue = scatter(ax, truth(par1), truth(par2), ...
                    100, o.TrueColor, '.');
                uistack(hTrue, 'top');
            end
        end

        set(ax, 'LineWidth',1.25, 'FontName','Arial', 'FontSize',6, ...
            'Layer','top', 'TickLength',[0 0]);
        box(ax, 'on');
    end
end
end

function samples = resamplePosterior(result, gene, ind, n)
%RESAMPLEPOSTERIOR Return 3-by-n log10 samples using the particle weights.
if ~isstruct(result) || ~isscalar(result) || ...
        ~isfield(result,'ks_t') || ~isfield(result,'ws_t')
    error('Each inference result must be a structure with ks_t and ws_t.');
end
if ndims(result.ks_t) > 4 || size(result.ks_t,4) ~= 3 || ...
        ndims(result.ws_t) > 3 || ...
        gene > size(result.ks_t,1) || ind > size(result.ks_t,2) || ...
        gene > size(result.ws_t,1) || ind > size(result.ws_t,2)
    error('Check the gene, iteration, and documented inference array shapes.');
end
np = size(result.ks_t,3);
particles = reshape(result.ks_t(gene,ind,:,:), np, 3);
weights = reshape(result.ws_t(gene,ind,:), [], 1);
validateattributes(particles, {'numeric'}, {'real','finite','positive'});
validateattributes(weights, {'numeric'}, ...
    {'real','finite','nonnegative','numel',np});
if ~any(weights > 0)
    error('The selected gene/iteration must have at least one positive weight.');
end
idx = randsample(np, n, true, weights);
samples = log10(particles(idx,:)).';
end
