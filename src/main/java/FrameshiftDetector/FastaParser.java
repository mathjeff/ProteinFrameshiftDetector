package FrameshiftDetector;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

public class FastaParser {
  public static Map<String, String> parseFasta(File file) throws IOException {
    String path = file.getPath();
    List<String> lines = Files.readAllLines(Paths.get(path));
    return parseFasta(lines, path);
  }

  public static Map<String, String> parseFasta(List<String> lines, String filename) {
    Map<String, String> entries = new HashMap<String, String>();

    lines.add(">unused"); // We only know that one entry is done when we see the next one, so we add an unused protein at the end

    String name = null;
    StringBuilder contentBuilder = null;
    for (String line: lines) {
      if (line.startsWith(">")) {
        if (name != null) {
          if (entries.containsKey(name)) {
            throw new RuntimeException("Name " + name + " found multiple times in " + filename);
          }

          entries.put(name, contentBuilder.toString());
        }
        int spaceIndex = line.indexOf(' ');
        if (spaceIndex > 0)
          name = line.substring(1, spaceIndex); // remove ">" and everything after the space
        else
          name = line.substring(1); // remove ">"
        contentBuilder = new StringBuilder();
      } else {
        contentBuilder.append(line.replace("\n", "").replace("*",""));
      }
    }
    return entries;
  }
}
