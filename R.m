function [ R ,Rnormed] = R( xn,xn1,prop)

[R(1:prop.NB,1), Rnormed(1:prop.NB,1)] = RmassFull(xn,xn1,prop);
[R(prop.NB+1:2*prop.NB,1), Rnormed(prop.NB+1:2*prop.NB,1)] = RenergyFull(xn,xn1,xn1,prop);


end

