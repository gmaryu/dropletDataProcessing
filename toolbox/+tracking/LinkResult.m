classdef LinkResult
    % TODO: Add description 
    properties
        assigned (:,2) {mustBeInteger, mustBeNonnegative}
        unassignedFromFirstInput (:,1) {mustBeInteger, mustBeNonnegative}
        unassignedFromSecondInput (:,1) {mustBeInteger, mustBeNonnegative}
    end
    methods
        function obj = LinkResult(assigned, unassignedFromFirstInput, unassignedFromSecondInput)
            % TODO: Add description
            obj.assigned = assigned;
            obj.unassignedFromFirstInput = unassignedFromFirstInput;
            obj.unassignedFromSecondInput = unassignedFromSecondInput;
        end
    end
end