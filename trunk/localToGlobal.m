% psi is the local coordinate within the element (-1 to 1)
% h is element width
% x0 is the global coordinate of node0 of this element
function x = localToGlobal(x0, psi, h)

x = x0 + (h*(1.0 + psi)/2.0);

