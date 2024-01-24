classdef LinkResult
    %   LinkResult Class for storing the information about the process of connecting two segmentation results together
    %
    %   Constructor:
    %       linkResult = tracking.LinkResult(assigned, unassignedFromFirstInput, unassignedFromSecondInput)
    properties
        assigned (:,2) {mustBeInteger, mustBeNonnegative}
        unassignedFromFirstInput (:,1) {mustBeInteger, mustBeNonnegative}
        unassignedFromSecondInput (:,1) {mustBeInteger, mustBeNonnegative}
    end
    methods
        function obj = LinkResult(assigned, unassignedFromFirstInput, unassignedFromSecondInput)
            %   Creates a LinkResult object
            %
            %   Inputs:
            %       assigned - Nx2 matrix of indices connecting labels in one segmentation result to labels in another segmentation result
            %       unassignedFromFirstInput - Nx1 vector of labels in the first segmentation result that were not connected to any labels in the second segmentation result
            %       unassignedFromSecondInput - Nx1 vector of labels in the second segmentation result that were not connected to any labels in the first segmentation result
            obj.assigned = assigned;
            obj.unassignedFromFirstInput = unassignedFromFirstInput;
            obj.unassignedFromSecondInput = unassignedFromSecondInput;
        end
    end
end