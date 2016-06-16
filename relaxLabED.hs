{-# LANGUAGE UnicodeSyntax #-}

import Data.Matrix
import Data.List
import Data.Maybe (catMaybes)

-- | Returns the list of coefficient weights
-- multiplied by the corresponding probability.
-- In this case, it acts as a neighbourhood function.
cᵢⱼ :: Eq a => (Int, Int) -> Matrix a -> [a]
cᵢⱼ (iˣ, iʸ) m 
  = catMaybes 
  $ map (\ (nˣ, nʸ) 
        -> safeGet nˣ nʸ m
        ) $ [(nˣ, nʸ) | nˣ <- [iˣ-1..iˣ+1], nʸ <- [iʸ-1..iʸ+1]] \\ [(iˣ, iʸ)]

-- | Returns all indices in a matrix.
allElems :: Matrix a -> [(Int,Int)]
allElems m = [(x, y) | x <- [1..ncols m], y <- [1..nrows m]]

-- | Applies the updating scheme for edge detection
-- to a given matrix and returns the result.
-- λ₁ denotes the 'edge' label and λ₂ denotes the 'not edge' label.
-- They also represent the corresponding edge probability matrices.
updatingScheme :: (Fractional a, Ord a, Num a) 
               => [Int] -> Matrix a -> Matrix a
updatingScheme r (λ₁ @ m₁)
  = fromList (nrows m₁) (ncols m₁)
  $ map (\ aᵢ
        ->   (pˢ aᵢ λ₁) * (qˢλ₁ aᵢ) 
        /  ( (pˢ aᵢ λ₁) * (qˢλ₁ aᵢ) + (pˢ aᵢ λ₂) * (qˢλ₂ aᵢ) ) 
        ) $ allElems m₁

  where pˢ aᵢ    = getElem (fst aᵢ) (snd aᵢ)
        qˢλ₁ aᵢ  = sum $  map (\ aⱼ -> aⱼ * fromIntegral (r !! 0)) (cᵢⱼ aᵢ m₁)
                       ++ map (\ aⱼ -> aⱼ * fromIntegral (r !! 1)) (cᵢⱼ aᵢ m₂)
        qˢλ₂ aᵢ  = sum $  map (\ aⱼ -> aⱼ * fromIntegral (r !! 2)) (cᵢⱼ aᵢ m₁)
                       ++ map (\ aⱼ -> aⱼ * fromIntegral (r !! 3)) (cᵢⱼ aᵢ m₂)
        λ₂ @ m₂  = fromList (nrows m₁) (ncols m₁) $ map (\ p -> 1 - p) (toList m₁)

-- | Applies the relaxation labeling to a matrix 'iterNum - 1' times
-- based on the 'r' compatibilities. Returns the result of the final iteration.
relaxation :: (Fractional a, Ord a, Num a) 
           => Matrix a -> [Int] -> Int -> Matrix a
relaxation m₁ r iterNum = last $ take iterNum $ iterate (updatingScheme r) m₁

main :: IO ()
main = do mapM_ (putStrLn . prettyMatrix) [ relaxation mᴬ [1,1,1,1] 3,
                                            relaxation mᴮ [2,1,1,1] 3,
                                            relaxation mᶜ [2,1,1,1] 3
                                          ]

mᴬ = fromLists [ [0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.1, 0.1, 0.1, 0.0],
                 [0.0, 0.1, 0.9, 0.0, 0.0],
                 [0.0, 0.1, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0]
               ]

mᴮ = fromLists [ [0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 1.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0]
               ]

mᶜ = fromLists [ [0.0, 0.0, 0.0, 0.0, 0.0],
                 [0.0, 0.1, 0.1, 0.1, 0.0],
                 [0.0, 0.1, 1.0, 0.1, 0.0],
                 [0.0, 0.1, 0.1, 0.1, 0.0],
                 [0.0, 0.0, 0.0, 0.0, 0.0]
               ]