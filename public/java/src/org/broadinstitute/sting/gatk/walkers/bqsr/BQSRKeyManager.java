package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BitSetUtils;

import java.util.*;

/**
 * This class provides all the functionality for the BitSet representation of the keys to the hash table of BQSR
 *
 * It also handles the event type "covariate" which is not exactly a covariate, but is added as a key to the hashmap. The Key Manager will
 * add the event type as a bitset to the end of the covariate bitset key. This way, it won't get int the way of masking the information
 * out of the key for the actual covariates, and having the covariates handle it. The key manager handles the event type.
 *
 * The keys represented by this key manager will always have the same order:
 *
 * RequiredCovariate1, RequiredCovariate2, ..., RequiredCovariateN, OptionalCovariate1, OptionalCovariateID, EventType
 * RequiredCovariate1, RequiredCovariate2, ..., RequiredCovariateN, OptionalCovariate2, OptionalCovariateID, EventType
 * ...
 * RequiredCovariate1, RequiredCovariate2, ..., RequiredCovariateN, OptionalCovariateN, OptionalCovariateID, EventType
 *
 *
 * Note that Optional Covariates are optional, and the Key Manager should operate without them if necessary.
 *
 * @author Mauricio Carneiro
 * @since 3/6/12
 */
public class BQSRKeyManager {

    private final List<Covariate> requiredCovariates;
    private final List<Covariate> optionalCovariates;
    private final List<RequiredCovariateInfo> requiredCovariatesInfo;
    private final List<OptionalCovariateInfo> optionalCovariatesInfo;
    private final Map<String, Short> covariateNameToIDMap;

    private int nRequiredBits;                                                                                          // Number of bits used to represent the required covariates
    private int nOptionalBits;                                                                                          // Number of bits used to represent the standard covaraites
    private final int nOptionalIDBits;                                                                                  // Number of bits used to represent the optional covariates IDs
    private final int totalNumberOfBits;                                                                                // Sum of all of the above plus the event bits
    
    private final BitSet optionalCovariateMask;                                                                         // Standard mask for optional covariates bitset
    private final BitSet optionalCovariateIDMask;                                                                       // Standard mask for optional covariates order bitset
    
    /**
     * Initializes the KeyManager with the total number of covariates to use
     *
     * @param requiredCovariates the ordered list of required covariates
     * @param optionalCovariates the ordered list of optional covariates
     */
    public BQSRKeyManager(List<Covariate> requiredCovariates, List<Covariate> optionalCovariates) {
        this.requiredCovariates = new ArrayList<Covariate>(requiredCovariates);
        this.optionalCovariates = new ArrayList<Covariate>(optionalCovariates);
        requiredCovariatesInfo = new ArrayList<RequiredCovariateInfo>(requiredCovariates.size());                       // initialize the required covariates list
        optionalCovariatesInfo = new ArrayList<OptionalCovariateInfo>(optionalCovariates.size());                       // initialize the optional covariates list (size may be 0, it's okay)
        covariateNameToIDMap = new HashMap<String, Short>(optionalCovariates.size()*2);                                 // the map from covariate name to covariate id (when reading GATK Reports, we get the IDs as names of covariates)
        
        nRequiredBits = 0;
        for (Covariate required : requiredCovariates) {                                                                 // create a list of required covariates with the extra information for key management
            int nBits = required.numberOfBits();                                                                        // number of bits used by this covariate
            BitSet mask = genericMask(nRequiredBits, nBits);                                                            // create a mask for this covariate
            requiredCovariatesInfo.add(new RequiredCovariateInfo(nRequiredBits, mask, required));                       // Create an object for this required covariate
            nRequiredBits += nBits;
        }

        short id = 0;
        nOptionalBits = 0;
        for (Covariate optional : optionalCovariates) {
            int nBits = optional.numberOfBits();                                                                        // number of bits used by this covariate
            nOptionalBits = Math.max(nOptionalBits, nBits);                                                             // optional covariates are represented by the number of bits needed by biggest covariate
            BitSet optionalID = BitSetUtils.bitSetFrom(id);                                                             // calculate the optional covariate ID for this covariate
            optionalCovariatesInfo.add(new OptionalCovariateInfo(optionalID, optional));                                // optional covariates have standardized mask and number of bits, so no need to store in the RequiredCovariateInfo object
            String covariateName = optional.getClass().getSimpleName().split("Covariate")[0];                           // get the name of the covariate (without the "covariate" part of it) so we can match with the GATKReport
            covariateNameToIDMap.put(covariateName, id);
            id++;
        }

        nOptionalIDBits = BitSetUtils.numberOfBitsToRepresent(optionalCovariates.size());                               // number of bits used to represent the covariate ID
        optionalCovariateMask = genericMask(nRequiredBits, nOptionalBits);                                              // the generic mask to extract optional covariate bits from the combined bitset
        optionalCovariateIDMask = genericMask(nRequiredBits + nOptionalBits, nOptionalIDBits);                          // the generic mask to extract optional covariate ID bits from the combined bitset
        totalNumberOfBits = nRequiredBits + nOptionalBits + nOptionalIDBits + bitsInEventType();                        // total number of bits used in the final key
    }

    /**
     * Generates one key per optional covariate.
     * 
     * Keys include all required covariates, the standard covariate and the event type.
     *
     * Example allKeys:
     * RG, QUAL, CYCLE, CONTEXT
     *
     * List of BitSets returned by this example (given eventType):
     * RG, QUAL, CYCLE, EVENT
     * RG, QUAL, CONTEXT, EVENT
     *
     * Note: If there are no optional covariates, only one bitset key will be returned with all the required covariates and the event type
     *
     * @param allKeys      The keys in bitset representation for each covariate
     * @param eventType The type of event described by this keyset (e.g. mismatches, insertions, deletions)
     * @return one key in bitset representation per covariate
     */
    public List<BitSet> bitSetsFromAllKeys(BitSet[] allKeys, EventType eventType) {
        List<BitSet> allBitSets = new LinkedList<BitSet>();                                                             // Generate one key per optional covariate

        BitSet eventBitSet = BitSetUtils.bitSetFrom(eventType.index);                                                   // create a bitset with the event type
        int eventTypeBitIndex = nRequiredBits + nOptionalBits + nOptionalIDBits;                                        // Location in the bit set to add the event type bits

        int covariateIndex = 0;
        BitSet requiredKey = new BitSet(nRequiredBits);                                                                 // This will be a bitset holding all the required keys, to replicate later on
        for (RequiredCovariateInfo infoRequired : requiredCovariatesInfo)
            addBitSetToKeyAtLocation(requiredKey, allKeys[covariateIndex++], infoRequired.bitsBefore);                  // Add all the required covariates to the key set

        for (OptionalCovariateInfo infoOptional : optionalCovariatesInfo) {
            BitSet covariateKey = allKeys[covariateIndex++];                                                            // get the bitset from all keys
            if (covariateKey == null)
                continue;                                                                                               // do not add nulls to the final set of keys.

            BitSet optionalKey = new BitSet(totalNumberOfBits);                                                         // create a new key for this optional covariate
            optionalKey.or(requiredKey);                                                                                // import all the required covariates
            addBitSetToKeyAtLocation(optionalKey, covariateKey, nRequiredBits);                                         // add the optional covariate right after the required covariates
            addBitSetToKeyAtLocation(optionalKey, infoOptional.covariateID, nRequiredBits + nOptionalBits);             // add the optional covariate ID right after the optional covarite
            addBitSetToKeyAtLocation(optionalKey, eventBitSet, eventTypeBitIndex);                                      // Add the event type
            allBitSets.add(optionalKey);                                                                                // add this key to the list of keys
        }

        if (optionalCovariatesInfo.size() == 0) {                                                                       // special case when we have no optional covariates, add the event type to the required key (our only key)
            addBitSetToKeyAtLocation(requiredKey, eventBitSet, eventTypeBitIndex);                                      // Add the event type
            allBitSets.add(requiredKey);                                                                                // add this key to the list of keys
        }

        return allBitSets;
    }

    /**
     * Generates one bitset key for the covariates represented in Object[] key
     *
     * The covariates will have the actual objects produced by the covariates (probably read from the recalibration data file)
     * and will contain all required covariates and one (or none) optional covariates. Therefore, the product is one bitset key, not many.
     *
     * Example key:
     * RG, QUAL, CYCLE, CYCLE_ID, EventType
     *
     * @param key list of objects produced by the required covariates followed by one or zero optional covariates.
     * @return a bitset key representing these objects. Bitset encryption is done using the covariate's interface.
     */
    public BitSet bitSetFromKey(Object[] key) {
        BitSet bitSetKey = new BitSet(totalNumberOfBits);
        
        int requiredCovariate = 0;
        for (RequiredCovariateInfo infoRequired : requiredCovariatesInfo) {
            BitSet covariateBitSet = infoRequired.covariate.bitSetFromKey(key[requiredCovariate++]);                    // create a bitset from the object key provided using the required covariate's interface
            addBitSetToKeyAtLocation(bitSetKey, covariateBitSet, infoRequired.bitsBefore);                              // add it to the bitset key
        }
        
        if (optionalCovariatesInfo.size() > 0) {
            int optionalCovariate = requiredCovariatesInfo.size();                                                      // the optional covariate index in the key array
            int covariateIDIndex = optionalCovariate + 1;                                                               // the optional covariate ID index is right after the optional covariate's
            int covariateID = parseCovariateID(key[covariateIDIndex]);                                                  // when reading the GATK Report the ID may come in a String instead of an index
            OptionalCovariateInfo infoOptional = optionalCovariatesInfo.get(covariateID);                               // so we can get the optional covariate information
            
            BitSet covariateBitSet = infoOptional.covariate.bitSetFromKey(key[optionalCovariate]);                      // convert the optional covariate key into a bitset using the covariate's interface
            addBitSetToKeyAtLocation(bitSetKey, covariateBitSet, nRequiredBits);                                        // add the optional covariate right after the required covariates
            addBitSetToKeyAtLocation(bitSetKey, infoOptional.covariateID, nRequiredBits + nOptionalBits);               // add the optional covariate ID right after the optional covarite
        }
        
        int eventIndex = key.length - 1;                                                                                // the event type is always the last key
        int eventTypeBitIndex = nRequiredBits + nOptionalBits + nOptionalIDBits;                                        // location in the bit set to add the event type bits
        BitSet eventBitSet = bitSetFromEvent((EventType) key[eventIndex]);                                              // get the bit set representation of the event type
        addBitSetToKeyAtLocation(bitSetKey, eventBitSet, eventTypeBitIndex);                                            // add the event type

        return bitSetKey;
    }

    /**
     * Covariate id can be either the covariate name (String) or the actual id (short). This method
     * finds it's type and converts accordingly to the short notation.
     *
     * @param id the string or short representation of the optional covariate id
     * @return the short representation of the optional covariate id.
     */
    private short parseCovariateID(Object id) {
        return (id instanceof String) ? covariateNameToIDMap.get(id.toString()) : (Short) id;
    }

    /**
     * Generates a key set of objects from a combined bitset key.
     *
     * Masks out each covariate independently and decodes their values (Object) into a keyset
     *
     * @param key the bitset representation of the keys
     * @return an object array with the values for each key
     */
    public List<Object> keySetFrom(BitSet key) {
        List<Object> objectKeys = new ArrayList<Object>();
        for (RequiredCovariateInfo info : requiredCovariatesInfo) {
            BitSet covariateBitSet = extractBitSetFromKey(key, info.mask, info.bitsBefore);                             // get the covariate's bitset
            objectKeys.add(info.covariate.keyFromBitSet(covariateBitSet));                                              // convert the bitset to object using covariate's interface
        }

        if (optionalCovariatesInfo.size() > 0) {
            BitSet covBitSet = extractBitSetFromKey(key, optionalCovariateMask, nRequiredBits);                         // mask out the covariate bit set
            BitSet idbs = extractBitSetFromKey(key, optionalCovariateIDMask, nRequiredBits + nOptionalBits);            // mask out the covariate order (to identify which covariate this is)
            short id = BitSetUtils.shortFrom(idbs);                                                                     // covert the id bitset into a short
            Covariate covariate = optionalCovariatesInfo.get(id).covariate;                                                 // get the corresponding optional covariate object
            objectKeys.add(covariate.keyFromBitSet(covBitSet));                                                         // add the optional covariate to the key set
            objectKeys.add(covariate.getClass().getSimpleName().split("Covariate")[0]);                                 // add the covariate name using the id
        }
        objectKeys.add(eventFromBitSet(key));                                                                           // add the event type object to the key set

        return objectKeys;
    }

    /**
     * Translates a masked bitset into a bitset starting at 0
     *
     * @return a list of the optional covariates
     */
    public List<Covariate> getRequiredCovariates() {
        return requiredCovariates;
    }

    public List<Covariate> getOptionalCovariates() {
        return optionalCovariates;
    }

    /**
     * Translates a masked bitset into a bitset starting at 0
     *
     * @param key the masked out bitset
     * @param n   the number of bits to chop
     * @return a translated bitset starting at 0 for the covariate machinery to decode
     */
    private BitSet chopNBitsFrom(BitSet key, int n) {
        BitSet choppedKey = new BitSet();
        for (int i = key.nextSetBit(0); i >= 0; i = key.nextSetBit(i + 1))
            choppedKey.set(i - n);                                                                                      // Set every bit translocated to the beginning of the BitSet
        return choppedKey;
    }

    /**
     * Creates a mask for the requested covariate to extract the relevant bitset from a combined bitset key
     *
     * @param leadingBits the index of the covariate in the ordered covariate list
     * @param nBits the number of bits needed by the Covariate to represent its values in BitSet form
     * @return the bitset relevant to the covariate
     */
    
    private BitSet genericMask(int leadingBits, int nBits) {
        BitSet mask = new BitSet(leadingBits + nBits);
        mask.set(leadingBits, leadingBits + nBits);
        return mask;
    }

    /**
     * Decodes the event type (enum) from the full bitset key
     *
     * @param fullKey the full key of all covariates + event type
     * @return the decoded event type.
     */
    private EventType eventFromBitSet(BitSet fullKey) {
        BitSet eventKey = new BitSet();
        int firstBitIndex = nRequiredBits + nOptionalBits + nOptionalIDBits;
        for (int i = fullKey.nextSetBit(firstBitIndex); i >= 0; i = fullKey.nextSetBit(i + 1))
            eventKey.set(i - firstBitIndex);
        return EventType.eventFrom(BitSetUtils.shortFrom(eventKey));
    }
    
    private BitSet bitSetFromEvent(EventType eventType) {
        return BitSetUtils.bitSetFrom(eventType.index);
    }

    private int bitsInEventType() {
        return BitSetUtils.numberOfBitsToRepresent(EventType.values().length);
    }

    private void addBitSetToKeyAtLocation(BitSet key, BitSet bitSet, int location) {
        for (int j = bitSet.nextSetBit(0); j >= 0; j = bitSet.nextSetBit(j + 1))
            key.set(j + location);                                                                                      // translate the bits set in the key to their corresponding position in the full key
    }
    
    private BitSet extractBitSetFromKey (BitSet key, BitSet mask, int leadingBits) {
        BitSet bitSet = (BitSet) key.clone();                                  
        bitSet.and(mask);
        return chopNBitsFrom(bitSet, leadingBits);
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof BQSRKeyManager))
            return false;

        BQSRKeyManager other = (BQSRKeyManager) o;
        if (this == other)
            return true;

        if (requiredCovariatesInfo.size() != other.requiredCovariatesInfo.size() ||
                optionalCovariatesInfo.size() != other.optionalCovariatesInfo.size())
            return false;

        for (int i = 0; i < requiredCovariates.size(); i++) {
            Covariate myRequiredCovariate = requiredCovariates.get(i);
            Covariate otherRequiredCovariate = other.requiredCovariates.get(i);
            String thisName = myRequiredCovariate.getClass().getSimpleName();
            String otherName = otherRequiredCovariate.getClass().getSimpleName();
            if (!thisName.equals(otherName))
                return false;
        }

        for (int i = 0; i < optionalCovariates.size(); i++) {
            Covariate myOptionalCovariate = optionalCovariates.get(i);
            Covariate otherOptionalCovariate = other.optionalCovariates.get(i);
            String thisName = myOptionalCovariate.getClass().getSimpleName();
            String otherName = otherOptionalCovariate.getClass().getSimpleName();
            if (!thisName.equals(otherName))
                return false;
        }

        return true;
    }


    /**
     * Aggregate information for each Covariate
     */
    class RequiredCovariateInfo {
        public final int bitsBefore;                                                                                    // number of bits before this covariate in the combined bitset key
        public final BitSet mask;                                                                                       // the mask to pull out this covariate from the combined bitset key ( a mask made from bitsBefore and nBits )
        public final Covariate covariate;                                                                               // this allows reverse lookup of the Covariates in order

        RequiredCovariateInfo(int bitsBefore, BitSet mask, Covariate covariate) {
            this.bitsBefore = bitsBefore;
            this.mask = mask;
            this.covariate = covariate;
        }
    }

    class OptionalCovariateInfo {
        public final BitSet covariateID;                                                                                // cache the covariate ID
        public final Covariate covariate;

        OptionalCovariateInfo(BitSet covariateID, Covariate covariate) {
            this.covariateID = covariateID;
            this.covariate = covariate;
        }
    }
    
}
