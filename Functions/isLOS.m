function idxLOSBoleanArray = isLOS(S, objectPosition, object2TUAVElevation, q_TUAV)

idxLOSBoleanArray =  ones(length(S), 1);
for buildingIndex = 1:length(S)

    ySum = sum(S(buildingIndex).cornerBuildings(:, 1));
    xPosition = objectPosition(1)>S(buildingIndex).cornerBuildings(1, 2);

    if xPosition
        % Builing at the left side
        if ySum > 0
            % building at positive half
            thetaRightCorner2Object = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(1, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2))); % y / x
            thetaLeftCorner2Object = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(2, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(2, 2)));
            thetaTUAV2Object = atan2((objectPosition(2) - q_TUAV(2)), (objectPosition(1) - q_TUAV(1)));
            if thetaTUAV2Object > thetaRightCorner2Object && thetaTUAV2Object < thetaLeftCorner2Object
            %     disp('TUAV is in the cone')
                % In order to know whether the TUAV and RIS has LoS, elevation angles
                % must be compared
                Object2Building = (objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2)) / cos(thetaTUAV2Object);
                Object2BuildingElevation = atan((S(buildingIndex).height-objectPosition(3))/Object2Building);
                idxLOSBoleanArray(buildingIndex) = object2TUAVElevation >= Object2BuildingElevation;
            end
            
        elseif ySum < 0
            % building at negative half
            thetaRightCorner2Object = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(1, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(2, 2))); % y / x
            thetaLeftCorner2Object = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(2, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2)));
            thetaTUAV2Object = atan2((objectPosition(2) - q_TUAV(2)), (objectPosition(1) - q_TUAV(1)));
            if thetaTUAV2Object > thetaRightCorner2Object && thetaTUAV2Object < thetaLeftCorner2Object
            %     disp('TUAV is in the cone')
                % In order to know whether the TUAV and RIS has LoS, elevation angles
                % must be compared
                Object2Building = (objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2)) / cos(thetaTUAV2Object);
                Object2BuildingElevation = atan((S(buildingIndex).height-objectPosition(3))/Object2Building);
                idxLOSBoleanArray(buildingIndex)  = object2TUAVElevation >= Object2BuildingElevation;
            end
    
        else
            % building at the center
            thetaRightCorner2Object = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(1, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2))); % y/x
            thetaLeftCorner2Object = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(2, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2)));
            thetaTUAV2Object = atan2((objectPosition(2) - q_TUAV(2)), (objectPosition(1) - q_TUAV(1)));
            if thetaTUAV2Object > thetaRightCorner2Object && thetaTUAV2Object < thetaLeftCorner2Object
                % cone shadow
                Object2Building = (objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2))  / cos(thetaTUAV2Object);
                Object2BuildingElevation = atan((S(buildingIndex).height-objectPosition(3))/Object2Building);
                idxLOSBoleanArray(buildingIndex) = object2TUAVElevation >= Object2BuildingElevation;
            end
    
        end

    else
        % right
        if ySum > 0
            thetaRightCorner2Object    = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(2, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2))); % y / x
            thetaLeftCorner2Object     = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(1, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(2, 2)));
            thetaTUAV2Object           = atan2((q_TUAV(2)-objectPosition(2)), (q_TUAV(1)-objectPosition(1)));
            if thetaTUAV2Object > thetaRightCorner2Object && thetaTUAV2Object < thetaLeftCorner2Object
                Object2Building = (S(buildingIndex).cornerBuildings(2, 2) - objectPosition(1)) / cos(thetaTUAV2Object);
                Object2BuildingElevation = atan((S(buildingIndex).height-objectPosition(3))/Object2Building);
                idxLOSBoleanArray(buildingIndex) = object2TUAVElevation >= Object2BuildingElevation;
            end

        elseif ySum < 0
            thetaRightCorner2Object    = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(2, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(2, 2))); % y / x
            thetaLeftCorner2Object     = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(1, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(1, 2)));
            thetaTUAV2Object           = atan2((q_TUAV(2)-objectPosition(2)), (q_TUAV(1)-objectPosition(1)));
            if thetaTUAV2Object > thetaRightCorner2Object && thetaTUAV2Object < thetaLeftCorner2Object
                Object2Building = (S(buildingIndex).cornerBuildings(2, 2) - objectPosition(1)) / cos(thetaTUAV2Object);
                Object2BuildingElevation = atan((S(buildingIndex).height-objectPosition(3))/Object2Building);
                idxLOSBoleanArray(buildingIndex) = object2TUAVElevation >= Object2BuildingElevation;
            end

        else
            thetaRightCorner2Object    = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(2, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(2, 2))); % y / x
            thetaLeftCorner2Object     = atan((objectPosition(2) - S(buildingIndex).cornerBuildings(1, 1))/(objectPosition(1) - S(buildingIndex).cornerBuildings(2, 2)));
            thetaTUAV2Object           = atan2((q_TUAV(2)-objectPosition(2)), (q_TUAV(1)-objectPosition(1)));
            if thetaTUAV2Object > thetaRightCorner2Object && thetaTUAV2Object < thetaLeftCorner2Object
                Object2Building = (S(buildingIndex).cornerBuildings(2, 2) - objectPosition(1)) / cos(thetaTUAV2Object);
                Object2BuildingElevation = atan((S(buildingIndex).height-objectPosition(3))/Object2Building);
                idxLOSBoleanArray(buildingIndex) = object2TUAVElevation >= Object2BuildingElevation;
            end

        end

    end
    
end
