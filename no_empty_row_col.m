% need to mod the path below !
g = genpath('cesar7f-BiMat-96c324b');
addpath(g);

rand('seed', 100);

% load '/Users/psy/Downloads/cesar7f-BiMat-96c324b/examples/data/phage_bacteria_matrices.mat'

% matrix = phage_bacteria_matrices.matrices{38};

matrix = csvread('bimat_input_matrix_no_empty_rowcol.csv');

bp = Bipartite(matrix);

pformat = PlotFormat();
pformat.use_labels = false;
pformat.back_color = [110,110,20]/255;
pformat.cell_color = 'white';
pformat.use_isocline_module = false;
bp.plotter.SetPlotFormat(pformat);

font_size = 16;

figure(); set(gcf,'Position', [38,64,1290,435]);
subplot(1,3,1); bp.plotter.PlotMatrix; title('Original','FontSize',font_size);
subplot(1,3,2); bp.plotter.PlotNestedMatrix; title('Nested','FontSize',font_size);xlabel('','FontSize',font_size+6);
subplot(1,3,3); bp.plotter.PlotModularMatrix; title('Modular','FontSize',font_size);


stest_modul = StatisticalTest.TEST_COMMUNITY_STRUCTURE(matrix,200,@NullModels.AVERAGE,@LeadingEigenvector);
stest_nest = StatisticalTest.TEST_NESTEDNESS(matrix,200,@NullModels.EQUIPROBABLE,@NestednessNODF);
