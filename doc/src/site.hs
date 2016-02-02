import System.Console.GetOpt
import System.Environment
import System.Exit
import System.IO.Error
import System.IO
import Control.Monad (liftM, foldM, ap, when, (>=>))
import Data.Maybe
import Data.List (last)
import Hakyll
import Text.Pandoc
import Text.Pandoc.Options
import Text.Pandoc.Walk (walkM)
import Text.CSL.Pandoc
import Data.Bits ((.|.))
import qualified Text.Regex.PCRE as RE
import qualified System.FilePath as P
import Debug.Trace

-- command line option parsing
-- (we want to specify source, destination directories from the command-line
-- instead of in the source itself, which is what Hakyll expects.)

putStrLn' s = putStr s >> when (null s || last s /= '\n') (putStr "\n")
hPutStrLn' t s = hPutStr t s >> when (null s || last s /= '\n') (hPutStr t "\n")

data Options = Options {
   optVerbose :: Bool,
   optCommand :: String, 
   optCache :: Maybe String,
   optDest :: Maybe String,
   optSource :: [String] }

parseArgs args = do
    when (length args < 1) $ ioError $ userError "missing command\n"
    let opts0 = defaultOptions { optCommand = (head args) }
    when (not $ elem (optCommand opts0) ["rebuild", "clean"]) $
        ioError $ userError "unrecognized command"

    case getOpt Permute options (tail args) of
        (actions,rest,[]) -> do
            when (length rest > 1) $ ioError $ userError "too many source directories"
            opts <- foldM (flip ($)) opts0 actions
            return $ opts { optSource = rest }
        (_,_,err:_) -> ioError $ userError err

    where
        defaultOptions = Options False "" Nothing Nothing []
        usageString = "usage: site COMMAND [OPTIONS] SRCDIR\n" ++
                      "COMMAND is one of: rebuild, clean\n"++
                      "Options\n" ++
                      "    -c,--cache=DIR    set Hakyll cache directory\n" ++
                      "    -d,--dest=DIR     set destination directory\n" ++
                      "    -h,--help         display this usage information\n" ++
                      "    -v,--verbose      detailed logging/messages\n"
        options = [
            Option "c" ["cache"]
                (ReqArg (\a o ->
                    if isNothing (optCache o)
                    then return $ o {optCache = Just a}
                    else ioError $ userError "--cache option must be given at most once")
                    "DIR")
                "specify Hakyll cache directory",
            Option "d" ["dest"]
                (ReqArg (\a o ->
                    if isNothing (optDest o)
                    then return $ o {optDest = Just a}
                    else ioError $ userError "--dest option must be given at most once")
                    "DIR")
                "specify destination directory",
            Option "h" ["help"]
                (NoArg (\o -> putStrLn' usageString >> exitWith ExitSuccess))
                "print usage information",
            Option "v" ["verbose"]
                (NoArg (\o -> return $ o { optVerbose = True }))
                "print usage information"
            ]

-- run Hakyll with build or clean options only, with
-- directories specifed via command line options

main = do
    opts <- catchIOError (getArgs >>= parseArgs) (\e -> do
        hPutStrLn' stderr $ "site: " ++ (ioeGetErrorString e)
        hPutStrLn  stderr "Try 'site --help' for usage information."
        exitWith $ ExitFailure 2)
        
    let dc = defaultConfiguration
    let config = dc {
        destinationDirectory = fromMaybe (destinationDirectory dc) (optDest opts),
        storeDirectory = fromMaybe (storeDirectory dc) (optCache opts),
        providerDirectory = fromMaybe (providerDirectory dc) (listToMaybe $ optSource opts) }
    
    let hakyllArgs = (optCommand opts):if optVerbose opts then ["-v"] else []
    withArgs hakyllArgs $ hakyllWith config rules

-- customize Pandoc compiler process:

-- 1. transform internal links to their corresponding target
regexUTF8 = RE.makeRegexOpts (RE.defaultCompOpt .|. RE.compUTF8) RE.defaultExecOpt

transformInternal :: Pandoc -> Compiler Pandoc
transformInternal = walkM rewriteLink
    where rewriteLink link@(Link a es (url,title)) = do
                         newUrl <- localTarget url
                         return $ Link a es (newUrl,title)
          rewriteLink x = return x
          isExternal = RE.match (regexUTF8 "^[A-Za-z][A-Za-z0-9+-.]*:")
          localTarget url | isExternal url = return url
                          | otherwise = do
                                target <- getRoute $ fromFilePath url
                                let absTarget = liftM (P.combine "/") target
                                return $ fromMaybe url absTarget

-- 2. Use processCites' to add citations described by in-file YAML block

transformCites :: Pandoc -> Compiler Pandoc
transformCites = unsafeCompiler . processCites'

-- 3. bypass Pandoc table column wrap sizes
pandocReaderOptions = defaultHakyllReaderOptions { readerColumns = 200 }
pandocWriterOptions = defaultHakyllWriterOptions { writerWrapText = WrapAuto }

pandocCompiler' = pandocCompilerWithTransformM pandocReaderOptions pandocWriterOptions
                        (transformInternal >=> transformCites)


-- rules

rules = do
    match (fromGlob "templates/*.html") $ do
        compile templateCompiler
    match (fromGlob "css/*.css") $ do
        route idRoute
        compile compressCssCompiler
    match (fromGlob "**.md") $ do
        route $ setExtension "html"
        compile $
            pandocCompiler'  >>=
            loadAndApplyTemplate (fromFilePath "templates/default.html") defaultContext >>=
            relativizeUrls

