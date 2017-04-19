/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.PSMDataStructure;

import SortedListLib.SortedList;
import java.util.Comparator;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class SortedProteinListMaxIniProb extends SortedList<ProtID> {

    public SortedProteinListMaxIniProb() {
        super(new Comparator<ProtID>() {
            @Override
            public int compare(ProtID x, ProtID y) {
                if (x.MaxIniProb == y.MaxIniProb) {
                    return Float.compare(x.Mass, y.Mass);
                }
                return Float.compare(y.MaxIniProb, x.MaxIniProb);
            }
        });
    }

}
