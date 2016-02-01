import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO.Error
import System.IO
import Control.Monad (liftM, foldM, ap, when)
import Data.Maybe
import Hakyll.Core.Configuration

data Options = Options {
   optCommand :: String, 
   optCache :: Maybe String,
   optDest :: Maybe String,
   optSource :: [String] }

parseArgs args = do
    when (length args < 1) $ ioError $ userError "missing command\n"
    let opts0 = defaultOptions { optCommand = (head args) }
    when (not $ elem (optCommand opts0) ["build", "clean"]) $
        ioError $ userError "unrecognized command"
--
    case getOpt Permute options (tail args) of
        (actions,rest,[]) -> do
            when (length rest > 1) $ ioError $ userError "too many source directories\n"
            opts <- foldM (flip ($)) opts0 actions
            return $ opts { optSource = rest }
        (_,_,err:_) -> ioError $ userError err

    where
        defaultOptions = Options "" Nothing Nothing []
        usageString = "usage: site COMMAND [OPTIONS] SRCDIR\n" ++
                      "COMMAND is one of: build, clean\n"++
                      "Options\n" ++
                      "    -c,--cache=DIR    set Hakyll cache directory\n" ++
                      "    -d,--dest=DIR     set destination directory\n" ++
                      "    -h,--help         display this usage information"
        options = [
            Option "c" ["cache"]
                (ReqArg (\a o ->
                    if isNothing (optCache o)
                    then return $ o {optCache = Just a}
                    else ioError $ userError "--cache option must be given at most once\n")
                    "DIR")
                "specify Hakyll cache directory",
            Option "d" ["dest"]
                (ReqArg (\a o ->
                    if isNothing (optDest o)
                    then return $ o {optDest = Just a}
                    else ioError $ userError "--dest option must be given at most once\n")
                    "DIR")
                "specify destination directory",
            Option "h" ["help"]
                (NoArg (\o -> putStrLn usageString >> exitWith ExitSuccess))
                "print usage information"
            ]

main = do
    opts <- catchIOError (getArgs >>= parseArgs) (\e -> do
        hPutStrLn stderr $ "site: " ++ (ioeGetErrorString e)
            ++ "Try 'site --help' for usage information."
        exitWith $ ExitFailure 2)
        
    let dc = defaultConfiguration
    let config = dc {
        destinationDirectory = fromMaybe (destinationDirectory dc) (optDest opts),
        storeDirectory = fromMaybe (storeDirectory dc) (optCache opts),
        providerDirectory = fromMaybe (providerDirectory dc) (listToMaybe $ optSource opts) }
    
    
