function [griddedCruiseDIC, griddedCruiseTMP, griddedCruiseSAL] = gridAnyCruiseDICandTemp(EXPOCODE, smoothingWindow, removeSparseCutoff, xVariable, xlon_lat, z_all, mask)

% This function calculates the DIC and Temperature data on a uniform 20m vertical grid, for any 
% cruise in the GLODAP dataset. It is similar to gridCruiseData, except it should work for any 
% hydrographic occupation, it will only return DIC and Temperature, and you feed it an expocodeno 
% instead of a year.
% As before, it will perform a vertical interpolation, smooth horizontally, and then horizontally
% interpolate. Will this work nicely everywhere? Probably not. I suspect I'm going to need to 
% rewrite a lot.
% xVariable should be either latitude or longitude: we will use it to choose how to grid the outputs.
% We need to join some cruises, so EXPOCODE may be a cell array of strings or a single string.



if ~strcmpi(xVariable,'longitude') & ~strcmpi(xVariable,'latitude')
		error("'xVariable' must be either 'longitude' or 'latitude' (case insensitive)");
end

[DIC, obsLon, obsZ] = extractCruiseObservations(EXPOCODE, xVariable, 'G2tco2');
[TMP, ~, ~] = extractCruiseObservations(EXPOCODE, xVariable, 'G2theta');
[SAL, ~, ~] = extractCruiseObservations(EXPOCODE, xVariable, 'G2salinity');



% Now we need to start thinking about the logical flow of what happened in gridCruiseData, so 
% that we can try to replicate it.
% Luckily, the algorithm is commented to hell. 
% What it does is create a grid of NaN's, of the specified depth range and number of longitudes.
% Then, we run through each depth box and longitude box, and find any cruise data within these
% boxes. The algorithm here seems a touch convoluted by whatever.
% We then clear out any columns which don't have above a cutoff number of values in them back 
% to NaN, and then vertically interpolate to fill in all the NaNs in the columns which have enough
% values in them. 
% We then remove all the blank columns again, smooth the data we do have,


griddedCruiseTMP = createInitialGrid(TMP,xlon_lat,z_all,obsLon,obsZ);
griddedCruiseDIC = createInitialGrid(DIC,xlon_lat,z_all,obsLon,obsZ);
griddedCruiseSAL = createInitialGrid(SAL,xlon_lat,z_all,obsLon,obsZ);

scatteredDIC = griddedCruiseDIC;
scatteredTMP = griddedCruiseTMP;
scatteredSAL = griddedCruiseSAL;


griddedCruiseDIC = fillVerticalNaNs(griddedCruiseDIC,xlon_lat,z_all);
griddedCruiseSAL = fillVerticalNaNs(griddedCruiseSAL,xlon_lat,z_all);
griddedCruiseTMP = fillVerticalNaNs(griddedCruiseTMP,xlon_lat,z_all);

griddedCruiseDIC = fillHorizontalNaNs(griddedCruiseDIC,smoothingWindow,xlon_lat,z_all);
griddedCruiseSAL = fillHorizontalNaNs(griddedCruiseSAL,smoothingWindow,xlon_lat,z_all);
griddedCruiseTMP = fillHorizontalNaNs(griddedCruiseTMP,smoothingWindow,xlon_lat,z_all);

griddedCruiseTMP(isnan(mask)) = NaN;
griddedCruiseDIC(isnan(mask)) = NaN;
griddedCruiseSAL(isnan(mask)) = NaN;

if removeSparseCutoff
    griddedCruiseDIC = reassimilateSparseColumns(griddedCruiseDIC,scatteredDIC,xlon_lat,z_all,mask,smoothingWindow);
    griddedCruiseTMP = reassimilateSparseColumns(griddedCruiseTMP,scatteredTMP,xlon_lat,z_all,mask,smoothingWindow);
    griddedCruiseSAL = reassimilateSparseColumns(griddedCruiseSAL,scatteredSAL,xlon_lat,z_all,mask,smoothingWindow);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Body Ends Here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [output] = removeSparseColumns(input,cutoffValue)

xlon = evalin('caller','xlon;');
numVertPoints = NaN(244,1);
output = input;

for i = 1:length(xlon)
    numVertPoints(i) = sum(~isnan(input(:,i)));
    if numVertPoints(i) <= cutoffValue
        output(:,i) = NaN;
    end
end



end

function [output] = reassimilateSparseColumns(input,scatteredGrid,xGrid,yGrid,mask,smoothingWindow)


% Merge the interpolated grid with the gridded scattered values put onto a
% grid as regularly as possible.
griddedCruiseData = input;
griddedCruiseData(isnan(griddedCruiseData) & ~isnan(scatteredGrid)) = scatteredGrid(isnan(griddedCruiseData) & ~isnan(scatteredGrid));

% Aaaaand repeat the previous process

% Fill in nans vertically
maxdepth = zeros(size(xGrid));
for i = 1:length(xGrid)
    vals = find(~isnan(griddedCruiseData(:,i)));
    if ~isempty(vals)
        if length(vals) >1
            b = interp1(yGrid(vals),griddedCruiseData(vals,i),yGrid(vals(1):vals(end)));
            griddedCruiseData(vals(1):vals(end),i) = b;  
            maxdepth(i) = yGrid(vals(end));
        end
        
        if vals(1) >1
            griddedCruiseData(1:vals(1)-1,i) = griddedCruiseData(vals(1),i);
        end
    end
    
end


gridClone = griddedCruiseData;
[~,idx] = find(any(~isnan(gridClone), 1));
gridClone(:, ~any(~isnan(gridClone), 1))=[];   
% Above 3 lines remove any blank columns, the idea being that smoothing
% before filling them will work better. Then NaN out smoothed values where 
% they we're initially NaN.
testSmoothedGrid = smoothdata(gridClone,2,'rloess',smoothingWindow);
griddedCruiseData(:,idx) = testSmoothedGrid;

% Repeat horizontally

for i = 1:length(yGrid)
    vals = find(~isnan(griddedCruiseData(i,:)));
    if ~isempty(vals)
        if length(vals) >1
            b = interp1(xGrid(vals),griddedCruiseData(i,vals),xGrid(vals(1):vals(end)));
            griddedCruiseData(i,vals(1):vals(end)) = b;  
        end
    end
end

griddedCruiseData(isnan(mask)) = NaN;

output = griddedCruiseData;

end

function [output] = createInitialGrid(variable, xGrid, yGrid, obsLon, obsZ)
		
output = NaN(length(yGrid),length(xGrid));
for i = 1:length(yGrid)-1 % Run down each water column
    for j = 1:length(xGrid)-1 % Run along each longitude
        x_idx = find(obsLon<xGrid(j+1) & obsLon>xGrid(j)); % Find indices within the column we're looking at
        y_idx = find(obsZ<yGrid(i+1) & obsZ>yGrid(i));   % Find indices within the depth range we're looking at
        
        % Make the x and y lists the same length by padding the end of the
        % shorter list with NaN's
        if length(x_idx) > length(y_idx)
            y_idx = vertcat(y_idx,NaN(length(x_idx) - length(y_idx),1));
        else 
            x_idx = vertcat(x_idx,NaN(length(y_idx) - length(x_idx),1));
        end
        
        a = horzcat(x_idx,y_idx); % These 3 lines find repeated values.
        [~,idx] = unique(a);      % These repeated values are the values in 
        a(idx) = [];              % the lon-depth box .
        
        if length(a) == 1         % Only one value? Fill the box with it.
            output(i,j) = variable(a);
        end
        if length(a) >1           %^More than one value? Take the mean in the box.
            output(i,j) = nanmean(variable(a));
        end
    end
end

end


function [output] = fillVerticalNaNs(griddedCruiseData, xGrid, yGrid)


maxdepth = zeros(size(xGrid));
for i = 1:length(xGrid)
    vals = find(~isnan(griddedCruiseData(:,i)));
    if ~isempty(vals)
        if length(vals) >1
            b = interp1(yGrid(vals),griddedCruiseData(vals,i),yGrid(vals(1):vals(end)));
            griddedCruiseData(vals(1):vals(end),i) = b;  
            maxdepth(i) = yGrid(vals(end));
        end
        
        if vals(1) >1
            griddedCruiseData(1:vals(1)-1,i) = griddedCruiseData(vals(1),i);
        end
    end
    
end
output = griddedCruiseData;
end


function [output] = fillHorizontalNaNs(griddedCruiseData,smoothingWindow,xGrid,yGrid)

gridClone = griddedCruiseData;
[~,idx] = find(any(~isnan(gridClone), 1));
gridClone(:, ~any(~isnan(gridClone), 1))=[];   
% Above 3 lines remove any blank columns, the idea being that smoothing
% before filling them will work better. Then NaN out smoothed values where 
% they we're initially NaN.
testSmoothedGrid = smoothdata(gridClone,2,'rloess',smoothingWindow);
testSmoothedGrid(isnan(gridClone)) = NaN;
griddedCruiseData(:,idx) = testSmoothedGrid;


for i = 1:length(yGrid)
    vals = find(~isnan(griddedCruiseData(i,:)));
    if ~isempty(vals)
        if length(vals) > 1
            b = interp1(xGrid(vals),griddedCruiseData(i,vals),xGrid(vals(1):vals(end)));
            griddedCruiseData(i,vals(1):vals(end)) = b;  
        end
    end
end

output = griddedCruiseData;
end

