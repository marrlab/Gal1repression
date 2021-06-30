function [p,tab,chi2,labels] = mediantest(varargin)
%
% [p,tab,chi2,labels] = mediantest(varargin) performes Mood's mediantest: Mood's Median Test compares the medians of
% two or more groups. The test counts how many observations in each group
% are greater than the global median for all groups together and
% calculuates Chi-square statistics on those obervations.
% Less powerful than Kruskal-Wallis test, but requires fewer assumptions.
%
%
% Information from [1]:
%
% Appropriate data for use with Mood's median test:
% 	- One-way data with two or more groups
%   - Dependent variable is ordinal, interval, or ratio
%   - Independent variable is a factor with levels indicating groups
%   - Observations between groups are independent.  That is, not paired or repeated measures data
%
% Hypotheses
%   - Null hypothesis:  The medians of values for each group are equal.
%   - Alternative hypothesis (two-sided): The medians of values for each group are not equal.
%
% Interpretation
% - Significant results can be reported as �There was a significant difference in the median values among groups.�
%
% [1] Mangiafico, S.S. 2016. Summary and Analysis of Extension Program Evaluation in R, version 1.15.0. rcompanion.org/handbook/. (Pdf version: rcompanion.org/documents/RHandbookProgramEvaluation.pdf.)
%
% Inputs can be entered in several ways:
%     - individual data groups [e.g. mediantest(dataset1,dataset2,dataset3)] in which each dataset is a vector of continous responses, or data vector and group
%     - 2-column matrix containing responses (data) and groups identifier
%     - table containing columns 'data' and 'group'
% TEST FOR TOOLBOX
V = ver;
VName = {V.Name};
hasStatToolBox = any(cell2mat(strfind(VName,'Statistics'))==1);
if hasStatToolBox
	inputData = varargin;
	
	if istable(inputData{1})
		t = inputData{1};
		data = t.data;
		groups = t.groups;
	elseif isnumeric(inputData{1})
		if numel(inputData)<2 % NOT ENOUGH GROUPS
			warning('Please specify at least two groups.'), return
		elseif numel(inputData)>2 % assume data entered as several input arguments
			% SORT DATA
			[data, groups] = sortData(inputData);
		elseif numel(inputData) == 2 % COULD BE TWO GROUPS OR DATA/GROUP
			isCat = cellfun(@iscategorical, inputData);
			if any(isCat) % IF ONE VECTOR IS CATEGORICAL ONE INPUT ARE DATA SECOND IS GROUP
				data = inputData{~isCat};
				groups = inputData{isCat};
			else
				[data, groups] = sortData(inputData);
			end
		end
	else
		warning('Data type not supported. Must be numeric or table.'), return
	end
	
	[tab,chi2,p,labels] = crosstab(data>median(data),groups);
	
else
	warning('This function requires the Matlab Statistics Toolbox.')
	return
end
end
function [data,groups] = sortData(inputData)
nGroups = numel(inputData);
g = cell(nGroups,1);
d = cell(nGroups,1);
for iG=1:nGroups
	g{iG} = repmat(iG,size(inputData{iG}(:)));
	d{iG} = inputData{iG}(:);
end
groups = cell2mat(g);
data = cell2mat(d);
end
