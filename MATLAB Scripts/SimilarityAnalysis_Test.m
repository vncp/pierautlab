    clear
close all
clc

resultspath = '/Users/vince/lab/pierautlab/OriginalTest/data';
cd(resultspath)

TestNameFile

Fieldnames = {'Treatment', 'AnimalID', 'NeuronID', 'Protocol'};

Binsize = [0.005, 0.010, 0.020, 0.040, 0.050, 0.060, 0.080, 0.100, 0.250, 0.500, 1, 2];

%Creates struct for SimilarityAnalysis.m and calls the function for each bin for each stim for each dataitem
for b = 1:length(Binsize)
PARAMS.binsize = Binsize(b);
	sttime1 = cputime;
	for catnum = 1:length(Files.Stim.File)
		PARAMS.stimfilespec = Files.Stim.File{catnum};
		for itemnum = 1:length(Files.Resp.File{catnum}(:))
			PARAMS.datafilespec = Files.Resp.File{catnum}{itemnum};
				sttime2 = cputime;
				disp('**************************')
				Cdata{b}{catnum, itemnum} = SimilarityAnalysis(PARAMS);
				Cdata{b}{catnum, itemnum}.RecFile = ['{' num2str(catnum) '}{' num2str(itemnum) '}'];
					if isfield(Files, Fieldnames)
                    Cdata{b}{catnum, itemnum}.Treatment = Files.Treatment;
                    Cdata{b}{catnum, itemnum}.AnimalID = Files.AnimalID{catnum}{itemnum};
                    Cdata{b}{catnum, itemnum}.NeuronID = Files.NeuronID{catnum}{itemnum};
                    Cdata{b}{catnum, itemnum}.Protocol = Files.Protocol{catnum};
                	end
                endtime2 = cputime;
                disp([catnum itemnum])
                disp(['Time to do correlations = ' num2str(endtime2-sttime2)] )
                close all
                disp('***********************')
        end
    end
end

cd(resultspath)
outname = ['SimilarityAnalysis_TestOutput']
save(outname)