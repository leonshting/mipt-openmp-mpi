
import System.Random
import Data.Fixed

freq = 4.0
lengthi = 1.0
numof = 100000

kX::(Double, Double, Double, Double) -> Double
kX (x1, x2, y1, y2) = (y2-y1)/(x2-x1)

closest :: Double -> Double -> Double
closest xc freq =
        if (xc `mod'` freq) > freq/2
        then freq - xc `mod'` freq
        else xc `mod'` freq

ifintersect :: (Double, Double, Double, Double) -> Double -> Bool
ifintersect (x1,x2, x3, x4) freq =
    if ((closest x1 freq)^2 + (kX (x1,x2,x3,x4) * (closest x1 freq))^2 ) > lengthi^2
    then False
    else True

to_iter :: Bool -> Int -> Int
to_iter b a =
    if b
    then a+1
    else a


process :: Int -> Int -> Double
process a b = 2*2*lengthi/freq/((fromIntegral a)/(fromIntegral b))

main = do

    putStrLn (show (process (foldr (\x y ->to_iter (ifintersect x (freq :: Double)) y) 0 (take numof [(i,j,k,l)| i<-randomRs ((-50.0, 50.0) :: (Double, Double)) (mkStdGen 12)| j<-randomRs ((-50.0, 50.0) :: (Double, Double)) (mkStdGen 43)| k<-randomRs ((-50.0, 50.0) :: (Double, Double)) (mkStdGen 45)| l<-randomRs ((-50.0, 50.0) :: (Double, Double)) (mkStdGen 47)])) numof))

